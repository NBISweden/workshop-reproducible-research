"use strict";

// The full page consists of a main window with navigation and table of contents, and an inner
// iframe containing the current article. Which article is shown is determined by the main
// window's #hash portion of the URL. In fact, we use the simple rule: main window's URL of
// "rootUrl#relPath" corresponds to iframe's URL of "rootUrl/relPath".

var outerWindow = is_outer_page ? window : (window === window.parent ? null : window.parent);
var innerWindow = null;
var rootUrl = qualifyUrl(base_url);
var searchIndex = null;
var showPageToc = true;
var MutationObserver = window.MutationObserver || window.WebKitMutationObserver;

var Keys = {
	ENTER:  13,
	ESCAPE: 27,
	UP:     38,
	DOWN:   40,
};
var TABLE_CLASSES = [ 'table', 'table-striped', 'table-hover', 'table-bordered', 'table-condensed' ];

function onReady(doc, fun) {
	if (doc.readyState === 'interactive' || doc.readyState === 'complete') {
		fun();
	} else {
		doc.addEventListener('DOMContentLoaded', fun);
	}
}

function forEach(iterable, fun) {
	for (var i = 0; i < iterable.length; i++) {
		fun(iterable[i]);
	}
}

function toggleCollapse(collapse, value) {
	if (value) {
		collapse.show();
	} else {
		collapse.hide();
	}
}

function getCollapse(ele) {
	var collapse = bootstrap.Collapse.getInstance(ele);
	return collapse ? collapse : new bootstrap.Collapse(ele, { toggle: false });
}

if (is_outer_page) {
	// Main window.
	onReady(document, () => {
		innerWindow = document.getElementsByClassName('wm-article')[0].contentWindow;
		initMainWindow();
		ensureIframeLoaded();
	});

} else {
	// Article contents.
	innerWindow = window;
	if (outerWindow) {
		outerWindow.onIframeLoad();
	} else {
		// This is a page that ought to be in an iframe. Redirect to load the top page instead.
		var topUrl = getAbsUrl('#', getRelPath('/', window.location.href));
		if (topUrl) {
			window.location.href = topUrl;
		}
	}

	// Other initialization of iframe contents.
	hljs.initHighlightingOnLoad();
	onReady(document, function() {
		forEach(document.getElementsByTagName('table'), table => {
			TABLE_CLASSES.forEach(function(clazz) { table.classList.add(clazz); });
		});
		forEach(document.querySelectorAll('a[href*="://"]'), link => {
			link.target = '_blank';
		});
	});
}

/**
 * Creates event handlers for click and the enter key
 */
function onActivate(sel, fun) {
	forEach(document.querySelectorAll(sel), ele => {
		ele.addEventListener('click', () => fun(ele));
		ele.addEventListener('keydown', e => {
			if (e.which === 13 || e.which === 32) {
				fun(ele);
			}
		});
	});
}

/**
 * Returns whether to use small-screen mode. Note that the same size is used in css @media block.
 */
function isSmallScreen() {
	return window.matchMedia("(max-width: 767.98px)").matches;
}

/**
 * Given a relative URL, returns the absolute one, relying on the browser to convert it.
 */
function qualifyUrl(url) {
	var a = document.createElement('a');
	a.href = url;
	return a.href;
}

/**
 * Turns an absolute path to relative, stripping out rootUrl + separator.
 */
function getRelPath(separator, absUrl) {
	var prefix = rootUrl + (rootUrl.endsWith(separator) ? '' : separator);
	return absUrl.startsWith(prefix) ? absUrl.slice(prefix.length) : null;
}

/**
 * Turns a relative path to absolute, adding a prefix of rootUrl + separator.
 */
function getAbsUrl(separator, relPath) {
	var sep = rootUrl.endsWith(separator) ? '' : separator;
	return relPath === null ? null : rootUrl + sep + relPath;
}

/**
 * Redirects the iframe to reflect the path represented by the main window's current URL.
 * (In our design, nothing should change iframe's src except via updateIframe(), or back/forward
 * history is likely to get messed up.)
 */
function updateIframe(enableForwardNav) {
	// Grey out the "forward" button if we don't expect 'forward' to work.
	var histBtn = document.getElementById('hist-fwd');
	if (histBtn) {
		histBtn.classList.toggle('bg-gray', !enableForwardNav);
	}

	var targetRelPath = (getRelPath('#', outerWindow.location.href) || "").replace("~", "#");
	var targetIframeUrl;
	if (targetRelPath.length <= 1) {
		targetIframeUrl = home_url;
	} else {
		targetIframeUrl = getAbsUrl('/', targetRelPath);
	}
	var innerLoc = innerWindow.location;
	var currentIframeUrl = _safeGetLocationHref(innerLoc);

	console.log("updateIframe: %s -> %s (%s)", currentIframeUrl, targetIframeUrl,
		currentIframeUrl === targetIframeUrl ? "same" : "replacing");

	if (currentIframeUrl !== targetIframeUrl) {
		innerLoc.replace(targetIframeUrl);
		onIframeBeforeLoad(targetIframeUrl);
	}
	document.body.scrollTop = 0;
}

/**
 * Returns location.href, catching exception that's triggered if the iframe is on a different domain.
 */
function _safeGetLocationHref(location) {
	try {
		return location.href;
	} catch (e) {
		return null;
	}
}

/**
 * Returns the value of the given parameter in the URL's query portion.
 */
function getParam(key) {
	return new URLSearchParams(document.location.search).get(key);
}

/**
 * Update the height of the iframe container. On small screens, we adjust it to fit the iframe
 * contents, so that the page scrolls as a whole rather than inside the iframe.
 */
function onResize() {
	if (isSmallScreen()) {
		var article = document.getElementsByClassName('wm-article')[0];
		if (article.getAttribute('scrolling') !== 'no') {
			document.getElementById('wm-search-form').classList.remove('show');
			article.setAttribute('scrolling', 'no');
		}
		document.getElementsByClassName('wm-content-pane')[0].style.height = innerWindow.document.body.offsetHeight + 20 + 'px';
	} else {
		document.getElementById('wm-search-form').classList.add('show');
		document.getElementsByClassName('wm-content-pane')[0].removeAttribute('style');
		document.getElementsByClassName('wm-article')[0].setAttribute('scrolling', 'auto');
	}
}

/**
 * Gets the dropdown for the search results menu
 */
function getSearchResultsDropdown() {
	return bootstrap.Dropdown.getOrCreateInstance(document.getElementById('mkdocs-search-results-toggle'));
}

/**
 * Close TOC on small screens and hide the search
 */
function closeTempItems() {
	if (isSmallScreen()) {
		forEach(document.getElementsByClassName('wm-toc-triggered'), ele => {
			ele.classList.remove('wm-toc-triggered');
		});
		getSearchResultsDropdown().hide();
	}
}

/**
 * Adjusts link to point to a top page, converting URL from "base/path#tag" to "base#path~tag".
 */
function adjustLink(linkEl) {
	var relPath = getRelPath('/', linkEl.href);
	if (relPath !== null) {
		var newUrl = getAbsUrl('#', relPath.replace('#', '~'));
		linkEl.href = newUrl;
	}
}

/**
 * Given a URL, strips query and fragment, returning just the path.
 */
function cleanUrlPath(relUrl) {
	return relUrl.replace(/[#?].*/, '');
}

/**
 * Initialize the main window.
 */
function initMainWindow() {
	// wm-toc-button either opens the table of contents in the side-pane, or (on smaller screens)
	// shows the side-pane as a drop-down.
	onActivate('#wm-toc-button', _ => {
		if (isSmallScreen()) {
			window.scroll(0,0);
		}
		document.getElementById('main-content').classList.toggle('wm-toc-triggered');
	});

	window.addEventListener('resize', onResize);
	window.addEventListener('blur', closeTempItems);

	// Connect up the Back and Forward buttons (if present).
	onActivate('#hist-back', _ => window.history.back());
	onActivate('#hist-fwd', _ => window.history.forward());

	// When we click on an opener in the table of contents, open it.
	onActivate('.wm-toc-pane .wm-toc-opener', ele => {
		ele.classList.toggle('wm-toc-open');
		getCollapse(ele.nextElementSibling).toggle();
	});
	onActivate('.wm-toc-li', ele => {
		if (!ele.classList.contains('wm-page-toc-opener')) {
			return;
		}
		// Ignore clicks while transitioning.
		if (ele.nextElementSibling.classList.contains('collapsing')) {
			return;
		}
		showPageToc = !showPageToc;
		ele.classList.toggle('wm-page-toc-open', showPageToc);
		toggleCollapse(getCollapse(ele.nextElementSibling), showPageToc);
	});
	onActivate('.wm-toc-pane a', () => {
		document.getElementsByClassName('wm-toc-pane')[0].classList.toggle('wm-toc-triggered');
	});

	// Once the article loads in the side-pane, close the dropdown.
	document.getElementsByClassName('wm-article')[0].addEventListener('load', () => {
		onInnerWindowUpdated();

		document.title = innerWindow.document.title;
		onResize();

		// We want to update content height whenever the height of the iframe's content changes.
		// Using MutationObserver seems to be the best way to do that.
		var observer = new MutationObserver(onResize);
		observer.observe(innerWindow.document.body, {
			attributes: true,
			childList: true,
			characterData: true,
			subtree: true
		});

		innerWindow.focus();
	});

	// Initialize search functionality.
	initSearch();

	// Load the iframe now, and whenever we navigate the top frame.
	setTimeout(() => updateIframe(false), 0);
	window.addEventListener('hashchange', () => updateIframe(false));
}

function onInnerWindowUpdated() {
	window.history.replaceState(null, '', getAbsUrl('#', getRelPath('/', innerWindow.location.href).replace('#', '~')));
	getSearchResultsDropdown().hide();
}

function onIframeBeforeLoad(url) {
	forEach(document.getElementsByClassName('wm-current'), ele => {
		ele.classList.remove('wm-current');
	});
	closeTempItems();

	var tocLi = getTocLi(url);
	if (tocLi) {
		tocLi.classList.add('wm-current');
	}
	for (; tocLi && !tocLi.classList.contains('wm-toc-pane'); tocLi = tocLi.parentElement) {
		if (tocLi.classList.contains('wm-toc-li-nested')) {
			tocLi.classList.remove('collapsing');
			tocLi.classList.add('collapse', 'show');
			tocLi.removeAttribute('style');
			tocLi.previousElementSibling.classList.add('wm-toc-open');
		}
	}
}

function getTocLi(url) {
	var relPath = getRelPath('/', cleanUrlPath(url));
	var elem = document.querySelector('.wm-article-link[href="' + relPath + '"]');
	if (elem) {
		return elem.parentElement;
	} else {
		return null;
	}
}

var _deferIframeLoad = false;

// Sometimes iframe is loaded before main window's ready callback. In this case, we defer
// onIframeLoad call until the main window has initialized.
function ensureIframeLoaded() {
	if (_deferIframeLoad) {
		onIframeLoad();
	}
}

function onIframeLoad() {
	if (!innerWindow) { _deferIframeLoad = true; return; }
	var url = innerWindow.location.href;
	onIframeBeforeLoad(url);
	innerWindow.addEventListener('hashchange', onInnerWindowUpdated);
	innerWindow.addEventListener('click', () => getSearchResultsDropdown().hide(), { capture: true });

	if (innerWindow.pageToc) {
		var relPath = getAbsUrl('#', getRelPath('/', cleanUrlPath(url)));
		var li = getTocLi(url);
		if (!li.classList.contains('wm-page-toc-open')) {
			renderPageToc(li, relPath, innerWindow.pageToc);
		}
	}

	innerWindow.focus();
}

/**
 * Hides a bootstrap collapsible element, and removes it from DOM once hidden.
 */
function collapseAndRemove(collapsibleElem) {
	collapsibleElem.addEventListener('hidden.bs.collapse', () => {
		collapsibleElem.parentElement.removeChild(collapsibleElem);
	});
	getCollapse(collapsibleElem).hide();
}

function renderPageToc(parentElem, pageUrl, pageToc) {
	var ul = document.createElement('ul');
	ul.classList.add('wm-page-toc-tree');
	function addItem(tocItem) {
		var li = document.createElement('li');
		li.classList.add('wm-toc-li');
		var link = document.createElement('a');
		link.href = pageUrl + tocItem.url.replace('#', '~');
		link.classList.add('wm-article-link', 'wm-page-toc-text');
		link.innerHTML = tocItem.title;
		li.append(link);
		ul.append(li);

		if (tocItem.children) {
			tocItem.children.forEach(addItem);
		}
	}
	pageToc.forEach(addItem);

	forEach(document.getElementsByClassName('wm-page-toc-opener'), ele => {
		ele.classList.remove('wm-page-toc-open');
		ele.classList.remove('wm-page-toc-opener');
	});
	forEach(document.getElementsByClassName('wm-page-toc'), collapseAndRemove);

	var li = document.createElement('li');
	li.classList.add('wm-page-toc', 'wm-toc-li-nested', 'collapse');
	li.append(ul);
	parentElem.insertAdjacentElement('afterend', li);

	parentElem.classList.add('wm-page-toc-opener');
	parentElem.classList.toggle('wm-page-toc-open', showPageToc);
	toggleCollapse(getCollapse(li), showPageToc);
}

var searchIndexReady = false;

/**
 * Initialize search functionality.
 */
function initSearch() {
	// Create elasticlunr index.
	searchIndex = elasticlunr(function() {
		this.setRef('location');
		this.addField('title');
		this.addField('text');
	});

	var searchBox = document.getElementById('mkdocs-search-query');
	var searchResults = document.getElementById('mkdocs-search-results');

	// Fetch the prebuilt index data, and add to the index.
	var req = new XMLHttpRequest();
	req.onreadystatechange = () => {
		if (req.readyState === 4) {
			if (req.status === 200) {
				var data = JSON.parse(req.responseText);
				data.docs.forEach(doc => {
					searchIndex.addDoc(doc);
				});
				searchIndexReady = true;
				document.dispatchEvent(new Event('searchIndexReady'));
			}
		}
	};
	req.open('GET', base_url + '/search/search_index.json', true);
	req.send();

	function showSearchResults(optShow) {
		var show = (optShow === false ? false : Boolean(searchBox.value));
		if (show) {
			doSearch({
				resultsElem: searchResults,
				query: searchBox.value,
				snippetLen: 100,
				limit: 10
			});
		}
		toggleCollapse(getSearchResultsDropdown(), show);
		return show;
	}

	searchBox.addEventListener('click', e => {
		if (!searchResults.classList.contains('show')) {
			if (showSearchResults()) {
				e.stopPropagation();
			}
		}
	});

	// Search automatically and show results on keyup event.
	searchBox.addEventListener('keyup', e => {
		var show = (e.which !== Keys.ESCAPE && e.which !== Keys.ENTER);
		showSearchResults(show);
	});

	// Open the search box (and run the search) on up/down arrow keys.
	searchBox.addEventListener('keydown', e => {
		if (e.which === Keys.UP || e.which === Keys.DOWN) {
			if (showSearchResults()) {
				e.stopPropagation();
				e.preventDefault();
				setTimeout(() => {
					var results = searchResults.getElementsByTagName('a');
					if (e.which === Keys.UP) {
						results[results.length - 1].focus();
					} else {
						results[0].focus();
					}
				}, 0);
			}
		}
	});

	searchResults.addEventListener('keydown', e => {
		if (e.which !== Keys.ENTER) {
			searchBox.focus();
		}
	});
}

function escapeRegex(s) {
	return s.replace(/[-\/\\^$*+?.()|[\]{}]/g, '\\$&');
}

/**
 * This helps construct useful snippets to show in search results, and highlight matches.
 */
function SnippetBuilder(query) {
	var termsPattern = elasticlunr.tokenizer(query).map(escapeRegex).join("|");
	this._termsRegex = termsPattern ? new RegExp(termsPattern, "gi") : null;
}

SnippetBuilder.prototype.getSnippet = function(text, len) {
	if (!this._termsRegex) {
		return text.slice(0, len);
	}

	// Find a position that includes something we searched for.
	var pos = text.search(this._termsRegex);
	if (pos < 0) { pos = 0; }

	// Find a period before that position (a good starting point).
	var start = text.lastIndexOf('.', pos) + 1;
	if (pos - start > 30) {
		// If too long to previous period, give it 30 characters, and find a space before that.
		start = text.lastIndexOf(' ', pos - 30) + 1;
	}
	var rawSnippet = text.slice(start, start + len);
	return rawSnippet.replace(this._termsRegex, '<b>$&</b>');
};

/**
 * Search the elasticlunr index for the given query, and populate the dropdown with results.
 */
function doSearch(options) {
	var resultsElem = options.resultsElem;
	resultsElem.innerHTML = '';

	// If the index isn't ready, wait for it, and search again when ready.
	if (!searchIndexReady) {
		resultsElem.innerHTML = '<a class="dropdown-item disabled">SEARCHING...</a>';
		document.addEventListener('searchIndexReady', () => doSearch(options), { once: true });
		return;
	}

	var query = options.query;
	var snippetLen = options.snippetLen;
	var limit = options.limit;

	if (query === '') { return; }

	var results = searchIndex.search(query, {
		fields: { title: { boost: 10 }, text: { boost: 1 } },
		expand: true,
		bool: "AND"
	});

	var snippetBuilder = new SnippetBuilder(query);
	if (results.length > 0){
		var len = Math.min(results.length, limit || Infinity);
		var lastDoc; // Used to check for duplicate results due to nesting of headings
		for (var i = 0; i < len; i++) {
			var doc = searchIndex.documentStore.getDoc(results[i].ref);

			if (lastDoc && lastDoc.location.startsWith(doc.location) && doc.text.includes(lastDoc.text)) {
				continue;
			}
			lastDoc = doc;

			var item = document.createElement('a');
			item.classList.add('dropdown-item');
			item.setAttribute('href', limit == 0 ? doc.location : pathJoin(base_url, doc.location));
			var header = document.createElement('h6');
			header.classList.add('dropdown-header');
			header.innerHTML = snippetBuilder.getSnippet(doc.title, snippetLen);
			var text = document.createElement('p');
			text.innerHTML = snippetBuilder.getSnippet(doc.text, snippetLen);
			item.append(header, text);
			resultsElem.append(item);
		}
		if (limit != 0) {
			forEach(resultsElem.getElementsByTagName('a'), adjustLink);
		}
		if (limit) {
			var divider = document.createElement('div');
			divider.classList.add('dropdown-divider');
			divider.setAttribute('role', 'separator');
			resultsElem.appendChild(divider);
			var allResults = document.createElement('a');
			allResults.classList.add('dropdown-item');
			allResults.id = 'search-show-all';
			allResults.setAttribute('target', 'article');
			allResults.href = base_url + '/search.html?q=' + query;
			allResults.innerHTML = 'SEE ALL RESULTS';
			resultsElem.append(allResults);
		}
	} else {
		resultsElem.innerHTML = '<span class="dropdown-item-text">NO RESULTS FOUND</span>';
	}
}

function pathJoin(prefix, suffix) {
	var nPrefix = prefix.endsWith('/') ? prefix.slice(0, -1) : prefix;
	var nSuffix = suffix.startsWith('/') ? suffix.slice(1) : suffix;
	return nPrefix + '/' + nSuffix;
}

A *markup language* is a system for annotating text documents in order to *e.g.*
define formatting. HTML, if you are familiar with that, is an example of a
markup language. HTML uses tags, such as:

```html
<h1> Heading </h1>
<h2> Sub-heading </h2>
<a href="www.webpage.com"> Link </a>
<ul>
  <li> List-item1 </li>
  <li> List-item2 </li>
  <li> List-item3 </li>
</ul>
```

*Markdown* is a lightweight markup language which uses plain-text syntax in
order to be as unobtrusive as possible, so that a human can easily read it. Look
at the following toy example:

```markdown
# Heading

A [link](http://example.com).

## Sub-heading

Text attributes _italic_, *italic*, **bold**, `monospace`.

### Another deeper heading

Bullet list:

  * apples
  * oranges
  * pears
```

A markdown document can be converted to other formats, such as HTML or PDF, for
viewing in a browser or a PDF reader; in fact, the page you are reading right
now is written in markdown. Markdown is somewhat ill-defined, and as a
consequence of that there exist many implementations and extensions. They share
most of the syntax, however, with various additions on top.

There are a lot more things you can do with markdown than what we show here.
Indeed, this entire course is mostly written in markdown! You can read more
about markdown [here](https://www.markdownguide.org/getting-started/).

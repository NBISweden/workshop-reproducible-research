// quarto revealjs-header extension

(function revealjsHeaderExtension() {
  "use strict";

  /**
   * Convert reveal lazy-image attributes to eager images for header logos.
   * Reveal can emit data-src; this normalizes it to src when needed.
   * @param {ParentNode} root
   */
  function hydrateLogoImages(root) {
    const logoImages = root.querySelectorAll(
      ".header-logo-left img, .header-logo-right img"
    );

    logoImages.forEach((img) => {
      const hasSrc = img.getAttribute("src");
      const dataSrc = img.getAttribute("data-src");

      if (!hasSrc && dataSrc) {
        img.setAttribute("src", dataSrc);
        img.removeAttribute("data-src");
      }
    });
  }

  /**
   * Place the generated header in the reveal container once.
   */
  function mountHeader() {
    const header = document.querySelector("div.reveal-header");
    const revealRoot = document.querySelector(".reveal");

    if (!header || !revealRoot) {
      return;
    }

    if (header.dataset.qrhMounted === "true") {
      return;
    }

    revealRoot.insertBefore(header, revealRoot.firstChild);
    header.dataset.qrhMounted = "true";
    hydrateLogoImages(header);
  }

  /**
   * Publish reveal's current zoom factor as a CSS custom property so the
   * logo CSS can size itself in lockstep with the (transform-scaled)
   * slide content instead of staying pinned to a constant real-world size.
   * @param {Reveal} reveal
   */
  function publishScale(reveal) {
    const scale = typeof reveal.getScale === "function" ? reveal.getScale() : 1;
    document.documentElement.style.setProperty(
      "--qrh-scale",
      String(scale || 1)
    );
  }

  /**
   * Initialize extension behavior after reveal and DOM are ready.
   */
  function initialize() {
    const reveal = window.Reveal;

    if (!reveal || typeof reveal.isReady !== "function") {
      return;
    }

    const onReady = () => {
      mountHeader();
      publishScale(reveal);

      if (typeof reveal.on === "function") {
        reveal.on("resize", () => publishScale(reveal));
      }
    };

    if (reveal.isReady()) {
      onReady();
      return;
    }

    if (typeof reveal.on === "function") {
      reveal.on("ready", onReady);
    }
  }

  window.addEventListener("load", initialize, { once: true });
})();

// Accordion collapse handler for revealjs (no Bootstrap JS dependency)
(function () {
  "use strict";

  function handleAccordionClick(e) {
    var button = e.currentTarget;
    var targetSelector = button.getAttribute("data-bs-target");
    if (!targetSelector) return;
    var target = document.querySelector(targetSelector);
    if (!target) return;

    var parentSelector = target.getAttribute("data-bs-parent");
    var parent = parentSelector ? document.querySelector(parentSelector) : null;
    var isOpening = !target.classList.contains("show");

    // Close siblings in the same accordion group
    if (parent && isOpening) {
      var openItems = parent.querySelectorAll(".accordion-collapse.show");
      for (var i = 0; i < openItems.length; i++) {
        if (openItems[i] !== target) {
          openItems[i].classList.remove("show");
          var siblingBtn = parent.querySelector(
            '[data-bs-target="#' + openItems[i].id + '"]'
          );
          if (siblingBtn) {
            siblingBtn.classList.add("collapsed");
            siblingBtn.setAttribute("aria-expanded", "false");
          }
        }
      }
    }

    // Toggle current item
    if (isOpening) {
      target.classList.add("show");
      button.classList.remove("collapsed");
      button.setAttribute("aria-expanded", "true");
    } else {
      target.classList.remove("show");
      button.classList.add("collapsed");
      button.setAttribute("aria-expanded", "false");
    }
  }

  function initAccordions() {
    var buttons = document.querySelectorAll(
      ".quarto-accordion .accordion-button"
    );
    for (var i = 0; i < buttons.length; i++) {
      buttons[i].addEventListener("click", handleAccordionClick);
    }
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", initAccordions);
  } else {
    initAccordions();
  }
})();

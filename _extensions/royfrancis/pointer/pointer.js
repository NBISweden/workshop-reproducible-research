var RevealPointer = (function () {
  "use strict";

  var DEFAULT_CONFIG = {
    key: "q",
    color: "red",
    pointerSize: 16,
    alwaysVisible: false,
    trail: false,
    trailDuration: 150,
    trailSampling: 2,
    trailMaxPoints: 80,
  };

  var KEY_CODE_BY_NAME = {
    backspace: 8,
    tab: 9,
    enter: 13,
    shift: 16,
    ctrl: 17,
    alt: 18,
    pausebreak: 19,
    capslock: 20,
    esc: 27,
    escape: 27,
    space: 32,
    pageup: 33,
    pagedown: 34,
    end: 35,
    home: 36,
    leftarrow: 37,
    uparrow: 38,
    rightarrow: 39,
    downarrow: 40,
    insert: 45,
    delete: 46,
  };

  var TRANSLATE_RE = /translate\(([-\d.]+)px,\s*([-\d.]+)px\)/;
  var SCALE_RE = /scale\(([-\d.]+)\)/;
  var MATRIX_RE = /matrix\(([^)]+)\)/;

  /**
   * Convert raw plugin config values into a validated configuration object.
   *
   * @param {object} rawPointerConfig Pointer config from Reveal/Quarto metadata.
   * @returns {object} Normalized config with defaults and constraints applied.
   */
  function normalizeConfig(rawPointerConfig) {
    var raw = rawPointerConfig || {};
    var key = typeof raw.key === "string" ? raw.key.toLowerCase().trim() : "";
    var pointerSize = Number(raw.pointerSize);
    var trailDuration = Number(raw.trailDuration);
    var trailSampling = Number(raw.trailSampling);
    var trailMaxPoints = Number(raw.trailMaxPoints);

    return {
      key: key || DEFAULT_CONFIG.key,
      color: typeof raw.color === "string" ? raw.color : DEFAULT_CONFIG.color,
      pointerSize:
        Number.isFinite(pointerSize) && pointerSize > 0
          ? pointerSize
          : DEFAULT_CONFIG.pointerSize,
      alwaysVisible:
        typeof raw.alwaysVisible === "boolean"
          ? raw.alwaysVisible
          : DEFAULT_CONFIG.alwaysVisible,
      trail: typeof raw.trail === "boolean" ? raw.trail : DEFAULT_CONFIG.trail,
      trailDuration:
        Number.isFinite(trailDuration) && trailDuration >= 0
          ? trailDuration
          : DEFAULT_CONFIG.trailDuration,
      trailSampling:
        Number.isFinite(trailSampling) && trailSampling >= 0
          ? trailSampling
          : DEFAULT_CONFIG.trailSampling,
      trailMaxPoints:
        Number.isFinite(trailMaxPoints) && trailMaxPoints >= 2
          ? Math.round(trailMaxPoints)
          : DEFAULT_CONFIG.trailMaxPoints,
    };
  }

  /**
   * Convert a key label into the key binding shape Reveal.js expects.
   *
   * @param {string} key Normalized key name.
   * @returns {{ key: string, keyCode: number|undefined }} Key binding fields.
   */
  function resolveKeyBinding(key) {
    var normalized = (key || "").toLowerCase();

    if (KEY_CODE_BY_NAME[normalized] != null) {
      return { key: normalized, keyCode: KEY_CODE_BY_NAME[normalized] };
    }

    if (normalized.length === 1) {
      var charCode = normalized.toUpperCase().charCodeAt(0);
      if (charCode >= 48 && charCode <= 90) {
        return { key: normalized, keyCode: charCode };
      }
    }

    return { key: DEFAULT_CONFIG.key, keyCode: DEFAULT_CONFIG.key.toUpperCase().charCodeAt(0) };
  }

  /**
   * Parse a CSS transform string into translate + scale components.
   * Supports the forms emitted by Reveal.js.
   *
   * @param {string} transformValue CSS transform value.
   * @returns {{ x: number, y: number, scale: number }} Parsed transform.
   */
  function parseTransform(transformValue) {
    var transform = transformValue || "";
    var translateMatch;
    var scaleMatch;
    var matrixMatch;
    var matrixValues;

    if (!transform || transform === "none") {
      return { x: 0, y: 0, scale: 1 };
    }

    translateMatch = TRANSLATE_RE.exec(transform);
    scaleMatch = SCALE_RE.exec(transform);
    if (translateMatch || scaleMatch) {
      return {
        x: translateMatch ? Number.parseFloat(translateMatch[1]) : 0,
        y: translateMatch ? Number.parseFloat(translateMatch[2]) : 0,
        scale: scaleMatch ? Number.parseFloat(scaleMatch[1]) : 1,
      };
    }

    matrixMatch = MATRIX_RE.exec(transform);
    if (matrixMatch) {
      matrixValues = matrixMatch[1].split(",").map(function (v) {
        return Number.parseFloat(v.trim());
      });
      if (matrixValues.length === 6 && Number.isFinite(matrixValues[0])) {
        return {
          x: Number.isFinite(matrixValues[4]) ? matrixValues[4] : 0,
          y: Number.isFinite(matrixValues[5]) ? matrixValues[5] : 0,
          scale: matrixValues[0] || 1,
        };
      }
    }

    return { x: 0, y: 0, scale: 1 };
  }

  /**
   * Decide if a new trail point should be appended given sampling distance.
   *
   * @param {{x:number,y:number}|null} lastPoint Last stored point.
   * @param {{x:number,y:number}} nextPoint Current pointer location.
   * @param {number} sampling Minimum point spacing in pixels.
   * @returns {boolean} Whether to append the point.
   */
  function shouldAppendTrailPoint(lastPoint, nextPoint, sampling) {
    var dx;
    var dy;

    if (!lastPoint) {
      return true;
    }

    dx = nextPoint.x - lastPoint.x;
    dy = nextPoint.y - lastPoint.y;

    return Math.sqrt(dx * dx + dy * dy) >= sampling;
  }

  /**
   * Evaluate a Catmull-Rom spline point for a 2D + time sample.
   *
   * @param {{x:number,y:number,time:number}} p0
   * @param {{x:number,y:number,time:number}} p1
   * @param {{x:number,y:number,time:number}} p2
   * @param {{x:number,y:number,time:number}} p3
   * @param {number} t Parametric position in [0, 1].
   * @returns {{x:number,y:number,time:number}} Interpolated point.
   */
  function catmullRomPoint(p0, p1, p2, p3, t) {
    var t2 = t * t;
    var t3 = t2 * t;

    return {
      x:
        0.5 *
        ((2 * p1.x) +
          (-p0.x + p2.x) * t +
          (2 * p0.x - 5 * p1.x + 4 * p2.x - p3.x) * t2 +
          (-p0.x + 3 * p1.x - 3 * p2.x + p3.x) * t3),
      y:
        0.5 *
        ((2 * p1.y) +
          (-p0.y + p2.y) * t +
          (2 * p0.y - 5 * p1.y + 4 * p2.y - p3.y) * t2 +
          (-p0.y + 3 * p1.y - 3 * p2.y + p3.y) * t3),
      time: p1.time + (p2.time - p1.time) * t,
    };
  }

  function createRevealPointerPlugin() {
    var config = DEFAULT_CONFIG;
    var keyBinding = resolveKeyBinding(DEFAULT_CONFIG.key);
    var pointerEnabled = false;
    var pointerEl = null;
    var trailCanvas = null;
    var trailCtx = null;
    var frameHandle = null;
    var trailPoints = [];
    var revealTransform = { x: 0, y: 0, scale: 1 };

    var pointerState = {
      x: 0,
      y: 0,
      clientX: 0,
      clientY: 0,
      isVisible: false,
      hasPosition: false,
    };

    /** @returns {string} The active Reveal transform value from body styles. */
    function currentTransformValue() {
      return window.getComputedStyle(document.body).transform || document.body.style.transform || "";
    }

    /** Recompute Reveal translation and scale so the pointer stays aligned with slides. */
    function updateRevealTransform() {
      revealTransform = parseTransform(currentTransformValue());
    }

    /** Update pointer element coordinates, visibility, and size. */
    function renderPointer() {
      var safeScale = revealTransform.scale || 1;
      var pointerScale = safeScale === 1 ? 1 : 1 / safeScale;

      if (!pointerEl) {
        return;
      }

      pointerEl.style.top = String((pointerState.y - revealTransform.y) / safeScale) + "px";
      pointerEl.style.left = String((pointerState.x - revealTransform.x) / safeScale) + "px";
      pointerEl.style.opacity = pointerState.isVisible ? "0.8" : "0";
      pointerEl.style.width = String(config.pointerSize * pointerScale) + "px";
      pointerEl.style.height = String(config.pointerSize * pointerScale) + "px";
    }

    /** Create the canvas used for pointer trails if it does not exist. */
    function ensureTrailCanvas() {
      if (trailCanvas) {
        return;
      }

      trailCanvas = document.createElement("canvas");
      trailCanvas.className = "pointer-trail";
      trailCanvas.style.opacity = pointerState.isVisible ? "1" : "0";
      document.body.appendChild(trailCanvas);
      trailCtx = trailCanvas.getContext("2d");
      resizeTrailCanvas();
    }

    /** Resize trail canvas to viewport dimensions with device pixel ratio handling. */
    function resizeTrailCanvas() {
      var dpr;

      if (!trailCanvas || !trailCtx) {
        return;
      }

      dpr = window.devicePixelRatio || 1;
      trailCanvas.width = Math.floor(window.innerWidth * dpr);
      trailCanvas.height = Math.floor(window.innerHeight * dpr);
      trailCanvas.style.width = String(window.innerWidth) + "px";
      trailCanvas.style.height = String(window.innerHeight) + "px";
      trailCtx.setTransform(dpr, 0, 0, dpr, 0, 0);
    }

    /** Remove trail points that are older than the configured duration. */
    function pruneTrail(now) {
      var cutoff = now - config.trailDuration;
      var i = 0;

      while (i < trailPoints.length && trailPoints[i].time < cutoff) {
        i += 1;
      }
      if (i > 0) {
        trailPoints.splice(0, i);
      }
    }

    /**
     * Build a spline-densified trail so curved pointer motion remains smooth.
     *
     * @param {Array<{x:number,y:number,time:number}>} points Input points.
     * @returns {Array<{x:number,y:number,time:number}>} Smoothed points.
     */
    function smoothedPoints(points) {
      var smooth = [];
      var i;
      var p0;
      var p1;
      var p2;
      var p3;
      var dx;
      var dy;
      var distance;
      var steps;
      var s;
      var t;

      if (!points.length) {
        return smooth;
      }

      if (points.length < 3) {
        return points.slice();
      }

      smooth.push(points[0]);
      for (i = 0; i < points.length - 1; i += 1) {
        p0 = i > 0 ? points[i - 1] : points[i];
        p1 = points[i];
        p2 = points[i + 1];
        p3 = i + 2 < points.length ? points[i + 2] : points[i + 1];

        dx = p2.x - p1.x;
        dy = p2.y - p1.y;
        distance = Math.sqrt(dx * dx + dy * dy);
        steps = Math.max(1, Math.min(12, Math.ceil(distance / 2)));

        for (s = 1; s <= steps; s += 1) {
          t = s / steps;
          smooth.push(catmullRomPoint(p0, p1, p2, p3, t));
        }
      }

      return smooth;
    }

    /** Draw the tapered pointer trail for the current animation frame. */
    function drawTrail(now) {
      var points;
      var pointCount;
      var maxWidth;
      var i;
      var pPrev;
      var pCurr;
      var pNext;
      var lastPoint;
      var age;
      var life;
      var halfWidth;
      var alpha;
      var dx;
      var dy;
      var length;
      var nx;
      var ny;
      var prevX;
      var prevY;
      var nextX;
      var nextY;
      var pointData;
      var prevData;
      var currData;
      var latestPoint;

      if (!config.trail || !trailCtx || !trailCanvas) {
        return;
      }

      trailCtx.clearRect(0, 0, window.innerWidth, window.innerHeight);

      if (!pointerState.isVisible || !pointerState.hasPosition || config.trailDuration === 0) {
        trailPoints = [];
        return;
      }

      pruneTrail(now);
      lastPoint = trailPoints.length ? trailPoints[trailPoints.length - 1] : null;
      latestPoint = {
        x: pointerState.clientX,
        y: pointerState.clientY,
        time: now,
      };
      if (shouldAppendTrailPoint(lastPoint, latestPoint, config.trailSampling)) {
        trailPoints.push(latestPoint);
        if (trailPoints.length > config.trailMaxPoints) {
          trailPoints.splice(0, trailPoints.length - config.trailMaxPoints);
        }
      }

      points = smoothedPoints(trailPoints);
      pointCount = points.length;
      if (pointCount < 2) {
        return;
      }

      maxWidth = Math.max(1, config.pointerSize * 0.75);
      pointData = [];

      for (i = 0; i < pointCount; i += 1) {
        pCurr = points[i];
        pPrev = i > 0 ? points[i - 1] : points[i];
        pNext = i < pointCount - 1 ? points[i + 1] : points[i];

        prevX = pCurr.x - pPrev.x;
        prevY = pCurr.y - pPrev.y;
        nextX = pNext.x - pCurr.x;
        nextY = pNext.y - pCurr.y;

        dx = prevX + nextX;
        dy = prevY + nextY;
        if (Math.abs(dx) < 0.001 && Math.abs(dy) < 0.001) {
          dx = pNext.x - pPrev.x;
          dy = pNext.y - pPrev.y;
        }

        length = Math.sqrt(dx * dx + dy * dy) || 1;
        nx = -dy / length;
        ny = dx / length;

        age = now - pCurr.time;
        life = Math.max(0, 1 - age / config.trailDuration);
        halfWidth = Math.max(0.2, maxWidth * life * 0.5);
        alpha = 0.65 * life;

        pointData.push({
          x: pCurr.x,
          y: pCurr.y,
          nx: nx,
          ny: ny,
          halfWidth: halfWidth,
          alpha: alpha,
        });
      }

      trailCtx.fillStyle = config.color;
      for (i = 1; i < pointData.length; i += 1) {
        prevData = pointData[i - 1];
        currData = pointData[i];

        trailCtx.globalAlpha = Math.min(prevData.alpha, currData.alpha);
        trailCtx.beginPath();
        trailCtx.moveTo(prevData.x + prevData.nx * prevData.halfWidth, prevData.y + prevData.ny * prevData.halfWidth);
        trailCtx.lineTo(currData.x + currData.nx * currData.halfWidth, currData.y + currData.ny * currData.halfWidth);
        trailCtx.lineTo(currData.x - currData.nx * currData.halfWidth, currData.y - currData.ny * currData.halfWidth);
        trailCtx.lineTo(prevData.x - prevData.nx * prevData.halfWidth, prevData.y - prevData.ny * prevData.halfWidth);
        trailCtx.closePath();
        trailCtx.fill();
      }
      trailCtx.globalAlpha = 1;
    }

    /** Animation loop for pointer and trail rendering. */
    function renderFrame(now) {
      renderPointer();
      drawTrail(now || performance.now());
      frameHandle = requestAnimationFrame(renderFrame);
    }

    /** Cancel active animation frame loop if running. */
    function stopFrame() {
      if (frameHandle != null) {
        cancelAnimationFrame(frameHandle);
        frameHandle = null;
      }
    }

    /** Handle pointer movement updates and recalculate Reveal transforms. */
    function onMouseMove(event) {
      pointerState.x = event.pageX;
      pointerState.y = event.pageY;
      pointerState.clientX = event.clientX;
      pointerState.clientY = event.clientY;
      pointerState.hasPosition = true;
      updateRevealTransform();
      renderPointer();
    }

    /** Keep trail canvas in sync with viewport size changes. */
    function onResize() {
      resizeTrailCanvas();
    }

    /**
     * Enable or disable pointer mode and related event listeners.
     *
     * @param {boolean} nextEnabled Target enabled state.
     */
    function setPointerEnabled(nextEnabled) {
      pointerEnabled = Boolean(nextEnabled);
      pointerState.isVisible = pointerEnabled;

      if (pointerEnabled) {
        document.addEventListener("mousemove", onMouseMove);
        window.addEventListener("resize", onResize);
        document.body.classList.add("no-cursor");
        if (config.trail) {
          ensureTrailCanvas();
        }
        if (trailCanvas) {
          trailCanvas.style.opacity = "1";
        }
        if (frameHandle == null) {
          frameHandle = requestAnimationFrame(renderFrame);
        }
      } else {
        document.removeEventListener("mousemove", onMouseMove);
        window.removeEventListener("resize", onResize);
        document.body.classList.remove("no-cursor");
        trailPoints = [];
        if (trailCanvas && trailCtx) {
          trailCtx.clearRect(0, 0, window.innerWidth, window.innerHeight);
          trailCanvas.style.opacity = "0";
        }
        stopFrame();
        renderPointer();
      }
    }

    /** Toggle pointer mode for keyboard activation workflows. */
    function togglePointer() {
      setPointerEnabled(!pointerEnabled);
    }

    /** Remove plugin DOM artifacts and listeners for safe reinitialization. */
    function cleanup() {
      setPointerEnabled(false);

      if (pointerEl && pointerEl.parentNode) {
        pointerEl.parentNode.removeChild(pointerEl);
      }
      pointerEl = null;

      if (trailCanvas && trailCanvas.parentNode) {
        trailCanvas.parentNode.removeChild(trailCanvas);
      }
      trailCanvas = null;
      trailCtx = null;
      trailPoints = [];
      pointerState.hasPosition = false;
    }

    return {
      id: "pointer",
      init: function (deck) {
        var pluginConfig = deck.getConfig() || {};

        config = normalizeConfig(pluginConfig.pointer);
        keyBinding = resolveKeyBinding(config.key);

        pointerEl = document.createElement("div");
        pointerEl.className = "cursor-dot";
        pointerEl.style.width = String(config.pointerSize) + "px";
        pointerEl.style.height = String(config.pointerSize) + "px";
        pointerEl.style.backgroundColor = config.color;
        if (config.alwaysVisible) {
          pointerEl.style.opacity = "0.8";
        }
        document.body.appendChild(pointerEl);

        if (config.alwaysVisible) {
          setPointerEnabled(true);
          return;
        }

        deck.addKeyBinding(
          {
            keyCode: keyBinding.keyCode,
            key: keyBinding.key,
          },
          function () {
            togglePointer();
          },
        );
      },
      destroy: function () {
        cleanup();
      },
    };
  }

  createRevealPointerPlugin.__internals = {
    normalizeConfig: normalizeConfig,
    resolveKeyBinding: resolveKeyBinding,
    parseTransform: parseTransform,
    shouldAppendTrailPoint: shouldAppendTrailPoint,
  };

  return createRevealPointerPlugin;
})();

if (typeof module !== "undefined" && module.exports) {
  module.exports = {
    normalizeConfig: RevealPointer.__internals.normalizeConfig,
    resolveKeyBinding: RevealPointer.__internals.resolveKeyBinding,
    parseTransform: RevealPointer.__internals.parseTransform,
    shouldAppendTrailPoint: RevealPointer.__internals.shouldAppendTrailPoint,
  };
}

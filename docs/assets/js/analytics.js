/** GoatCounter helpers. No-op until a site code is set in index/gene/about. */

export function track(path, title) {
  if (window.goatcounter && typeof window.goatcounter.count === "function") {
    window.goatcounter.count({ path, title: title || path, event: true });
  }
}

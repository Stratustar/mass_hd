// Lazy attach: 100+ <video> elements with real sources would have the browser fetching a
// quarter of a gigabyte on load. src is set only when a card comes near the viewport.
const io = new IntersectionObserver((entries) => {
  for (const e of entries) {
    if (!e.isIntersecting) continue;
    const v = e.target;
    if (!v.src && v.dataset.src) { v.src = v.dataset.src; v.preload = "metadata"; }
    io.unobserve(v);
  }
}, { rootMargin: "600px 0px" });
document.querySelectorAll("video[data-src]").forEach(v => io.observe(v));

// Row sync: play the four starts at one tau_m together, from a common t = 0. Toggling off
// pauses and rewinds so the next press starts them level again.
document.querySelectorAll("[data-row]").forEach(row => {
  const btn = row.querySelector(".sync");
  if (!btn) return;
  const vids = () => [...row.querySelectorAll("video")];
  btn.addEventListener("click", () => {
    const on = btn.hasAttribute("data-on");
    if (on) {
      btn.removeAttribute("data-on");
      btn.textContent = "▶ 同步播放这一行";
      vids().forEach(v => { v.pause(); v.currentTime = 0; });
      return;
    }
    btn.setAttribute("data-on", "");
    btn.textContent = "■ 停止";
    vids().forEach(v => {
      if (!v.src && v.dataset.src) v.src = v.dataset.src;
      v.currentTime = 0;
      v.play().catch(() => {});
    });
  });
});

// Click to zoom. The dialog gets its own element rather than moving the card's, so the
// grid keeps its layout while the overlay is open.
const dlg = document.getElementById("zoom");
const dv = dlg.querySelector("video");
const dt = dlg.querySelector(".zt");
document.querySelectorAll(".card video").forEach(v => {
  v.addEventListener("click", () => {
    const fig = v.closest("figure");
    dv.src = v.src || v.dataset.src;
    dt.textContent = fig.querySelector("figcaption").textContent.replace(/\s+/g, " ").trim();
    dlg.showModal();
    dv.play().catch(() => {});
  });
});
dlg.addEventListener("close", () => { dv.pause(); dv.removeAttribute("src"); dv.load(); });
dlg.addEventListener("click", (e) => { if (e.target === dlg) dlg.close(); });

/*
    kick-geometry.js

    Interactive kick-geometry explorer for the kick_info section of the
    output docs. It draws (1) a top view of the orbit before and after the
    supernova and (2) a rotatable 3D view of how the supernova reorients the
    orbit, with the Euler angles that COSMIC stores in kick_info. The physics
    lives in kick-physics.js, which must be loaded first.
*/
(function () {
    "use strict";

    const P = window.CosmicKickPhysics;
    const D = Math.PI / 180;
    // V_c = (G M / a)^(1/2) in km/s, for M in Msun and a in Rsun
    const VC_SUN = Math.sqrt(1.32712440018e11 / 695700);

    const wrap180 = a => ((a + 180) % 360 + 360) % 360 - 180;
    const wrap360 = a => ((a % 360) + 360) % 360;
    const fmtAng = a => `${a < 0 ? "−" : ""}${Math.abs(a).toFixed(1)}°`;
    const SLIDERS = ["ecc", "M", "loga", "theta", "phi", "vk", "m1", "m2", "dm"];

    function initExplorer(root) {
        const q = sel => root.querySelector(sel);
        const el = Object.fromEntries(SLIDERS.map(id => [id, q(`[data-kg="${id}"]`)]));
        const out = id => q(`[data-kg-out="${id}"]`);
        const set = (id, text) => { const node = out(id); if (node) node.textContent = text; };
        const topLayer = q('[data-kg-layer="top"]');
        const svg3d = q('[data-kg-svg="3d"]');
        const layer3d = q('[data-kg-layer="3d"]');
        const showSecond = q('[data-kg="second"]');

        const DEFAULT_VIEW = { az: -58 * D, el: 24 * D };
        const view = { ...DEFAULT_VIEW };

        function state() {
            const m1 = +el.m1.value, m2 = +el.m2.value, dm = Math.min(+el.dm.value, +el.dm.max);
            const Mi = m1 + m2, a = Math.pow(10, +el.loga.value);
            const vc = VC_SUN * Math.sqrt(Mi / a);
            return {
                e: +el.ecc.value, M: +el.M.value * D, th: +el.theta.value * D, ph: +el.phi.value * D,
                m1, m2, dm, m1n: m1 - dm, Mi, mf: (Mi - dm) / Mi,
                a, vc, vkkm: +el.vk.value, vk: +el.vk.value / vc,
            };
        }

        function compute() {
            const s = state();
            const pre = P.preSupernova(s.e, s.M);
            const kick = P.kickVector(s.vk, s.th, s.ph);
            const post = P.postSupernova(pre, kick, s.mf);
            const euler = P.eulerAngles(pre, post);
            const vel = P.velocities(pre, post, kick, s.m1, s.m2, s.m1n);
            return { s, pre, kick, post, euler, vel };
        }

        // ---- top view (orbital plane seen from +z) -------------------------

        const SCALE = 100, OX = 330, OY = 205, VSCALE = 1;
        const sx = x => (OX + SCALE * x).toFixed(1);
        const sy = y => (OY - SCALE * y).toFixed(1);

        // counter-clockwise arc (in world coordinates) about (cx, cy)
        function arc(cx, cy, rad, a0, a1) {
            const d = ((a1 - a0) % (2 * Math.PI) + 2 * Math.PI) % (2 * Math.PI);
            if (d < 1e-3) return "";
            const x0 = cx + rad * Math.cos(a0), y0 = cy + rad * Math.sin(a0);
            const x1 = cx + rad * Math.cos(a0 + d), y1 = cy + rad * Math.sin(a0 + d);
            return `M ${sx(x0)} ${sy(y0)} A ${SCALE * rad} ${SCALE * rad} 0 ${d > Math.PI ? 1 : 0} 0 ${sx(x1)} ${sy(y1)}`;
        }

        // points along the post-SN conic (bound or not), in 3D
        function conicPoints(post, n, rmax) {
            const pts = [];
            if (post.hmag < 1e-9) return pts;
            const p = post.hmag * post.hmag / post.mu;
            const eh = post.ecc > 1e-6 ? P.unit(post.evec) : P.unit(cross0(post.h));
            const qh = P.cross(P.unit(post.h), eh);
            const fmax = post.ecc < 1 ? Math.PI : Math.acos(-1 / post.ecc) - 1e-3;
            for (let i = 0; i <= n; i++) {
                const f = -fmax + (2 * fmax * i) / n;
                const rr = p / (1 + post.ecc * Math.cos(f));
                if (!(rr > 0) || rr > rmax) { pts.push(null); continue; }
                pts.push(P.add(P.scale(eh, rr * Math.cos(f)), P.scale(qh, rr * Math.sin(f))));
            }
            return pts;
        }
        // any vector perpendicular to h (used when the new orbit is circular)
        const cross0 = h => (Math.abs(h[0]) < 0.9 ? P.cross(h, [1, 0, 0]) : P.cross(h, [0, 1, 0]));

        function renderTop(c) {
            const { s, pre, kick, post } = c;
            let svg = "";
            // axes
            svg += `<line x1="0" y1="${sy(0)}" x2="560" y2="${sy(0)}" style="stroke: var(--kg-rule)" stroke-dasharray="4 4"/>`;
            svg += `<text x="552" y="${OY - 8}" text-anchor="end" class="small">+x, periapsis direction</text>`;
            svg += `<g transform="translate(26 378)"><line x1="0" y1="0" x2="34" y2="0" style="stroke: var(--kg-muted)" stroke-width="1.2" marker-end="url(#kg-ah-ink)"/><line x1="0" y1="0" x2="0" y2="-34" style="stroke: var(--kg-muted)" stroke-width="1.2" marker-end="url(#kg-ah-ink)"/><text x="40" y="4" class="small">x</text><text x="-4" y="-40" class="small">y</text><circle cx="0" cy="0" r="4" fill="none" style="stroke: var(--kg-muted)"/><circle cx="0" cy="0" r="1.2" style="fill: var(--kg-muted)"/><text x="10" y="18" class="small">z, h out of page</text></g>`;

            // orbit before
            let d = "";
            for (let i = 0; i <= 240; i++) {
                const E = (i / 240) * 2 * Math.PI;
                d += `${i ? "L" : "M"} ${sx(Math.cos(E) - s.e)} ${sy(pre.q * Math.sin(E))} `;
            }
            svg += `<path d="${d}Z" fill="none" style="stroke: var(--kg-ink)" stroke-width="1.6"/>`;
            svg += `<circle cx="${sx(1 - s.e)}" cy="${sy(0)}" r="3" style="fill: var(--kg-ink)"/>`;

            // orbit after (projected onto the original plane)
            d = "";
            let pen = false;
            for (const pt of conicPoints(post, 720, 40)) {
                if (!pt) { pen = false; continue; }
                d += `${pen ? "L" : "M"} ${sx(pt[0])} ${sy(pt[1])} `;
                pen = true;
            }
            svg += `<path d="${d}" fill="none" clip-path="url(#kg-clip-top)" style="stroke: var(--kg-post)" stroke-width="1.6" stroke-dasharray="6 4"/>`;
            if (post.ecc > 0.02 && post.hmag > 1e-9) {
                const rp = post.hmag * post.hmag / post.mu / (1 + post.ecc);
                const eh = P.unit(post.evec);
                const px = rp * eh[0], py = rp * eh[1];
                svg += `<circle cx="${sx(px)}" cy="${sy(py)}" r="4" fill="none" style="stroke: var(--kg-post)" stroke-width="1.6"/>`;
                svg += `<text x="${(+sx(px) + 8).toFixed(1)}" y="${(+sy(py) + 16).toFixed(1)}" class="post" style="font-size: 11.5px">new periapsis</text>`;
            }

            // companion, separation and true anomaly
            const [rx, ry] = pre.r;
            svg += `<line x1="${sx(0)}" y1="${sy(0)}" x2="${sx(rx)}" y2="${sy(ry)}" style="stroke: var(--kg-muted)" stroke-width="1.2"/>`;
            svg += `<path d="${arc(0, 0, 0.3, 0, pre.nu)}" fill="none" style="stroke: var(--kg-muted)" stroke-width="1.1"/>`;
            const nl = pre.nu / 2 + (pre.nu < 0 ? Math.PI : 0);
            svg += `<text x="${sx(0.42 * Math.cos(nl))}" y="${(+sy(0.42 * Math.sin(nl)) + 4).toFixed(1)}" text-anchor="middle" class="small">ν</text>`;
            svg += `<circle cx="${sx(0)}" cy="${sy(0)}" r="7" style="fill: var(--kg-ink)"/>`;
            svg += `<text x="${sx(0)}" y="${OY + 24}" text-anchor="middle" class="small">companion</text>`;

            // theta, measured from +x at the star
            svg += `<line x1="${sx(rx)}" y1="${sy(ry)}" x2="${sx(rx + 0.5)}" y2="${sy(ry)}" style="stroke: var(--kg-muted)" stroke-width="1" stroke-dasharray="3 3"/>`;
            svg += `<path d="${arc(rx, ry, 0.26, 0, s.th)}" fill="none" style="stroke: var(--kg-kick)" stroke-width="1.1"/>`;
            svg += `<text x="${sx(rx + 0.38 * Math.cos(s.th / 2))}" y="${(+sy(ry + 0.38 * Math.sin(s.th / 2)) + 4).toFixed(1)}" text-anchor="middle" class="kick" style="font-size: 12px">θ</text>`;

            // orbital velocity
            const vx = rx + VSCALE * pre.v[0], vy = ry + VSCALE * pre.v[1];
            svg += `<line x1="${sx(rx)}" y1="${sy(ry)}" x2="${sx(vx)}" y2="${sy(vy)}" style="stroke: var(--kg-vel)" stroke-width="2.4" marker-end="url(#kg-ah-vel)"/>`;
            const vl = Math.hypot(pre.v[0], pre.v[1]);
            svg += `<text x="${sx(vx + 0.14 * pre.v[0] / vl)}" y="${(+sy(vy + 0.14 * pre.v[1] / vl) + 4).toFixed(1)}" text-anchor="middle" class="vel">v</text>`;

            // in-plane part of the kick
            const kin = Math.hypot(kick[0], kick[1]);
            if (kin > 1e-3) {
                const ks = Math.min(VSCALE, 2.0 / kin);
                const kx = rx + ks * kick[0], ky = ry + ks * kick[1];
                svg += `<line x1="${sx(rx)}" y1="${sy(ry)}" x2="${sx(kx)}" y2="${sy(ky)}" style="stroke: var(--kg-kick)" stroke-width="2.4" marker-end="url(#kg-ah-kick)"/>`;
                svg += `<text x="${sx(kx + 0.2 * kick[0] / kin)}" y="${(+sy(ky + 0.2 * kick[1] / kin) + 4).toFixed(1)}" text-anchor="middle" class="kick">${ks < VSCALE ? "kick (shortened)" : "kick"}</text>`;
            }
            svg += `<circle cx="${sx(rx)}" cy="${sy(ry)}" r="6" style="fill: var(--kg-panel); stroke: var(--kg-ink)" stroke-width="2"/>`;
            topLayer.innerHTML = svg;
        }

        // ---- 3D view (orthographic, rotatable) -----------------------------

        const S3 = 112, CX = 280, CY = 205;
        function projector() {
            const ca = Math.cos(view.az), sa = Math.sin(view.az);
            const ce = Math.cos(view.el), se = Math.sin(view.el);
            const right = [-sa, ca, 0], up = [-se * ca, -se * sa, ce];
            return p => [CX + S3 * P.dot(p, right), CY - S3 * P.dot(p, up)];
        }

        function render3D(c) {
            const { pre, post, euler } = c;
            const proj = projector();
            const pt = p => proj(p).map(v => v.toFixed(1)).join(" ");
            const line = (a, b, style, extra = "") =>
                `<line x1="${proj(a)[0].toFixed(1)}" y1="${proj(a)[1].toFixed(1)}" x2="${proj(b)[0].toFixed(1)}" y2="${proj(b)[1].toFixed(1)}" style="${style}" ${extra}/>`;
            const path = (pts, style, extra = "") => {
                let d = "", pen = false;
                for (const p of pts) {
                    if (!p) { pen = false; continue; }
                    d += `${pen ? "L" : "M"} ${pt(p)} `;
                    pen = true;
                }
                return `<path d="${d}" style="${style}" ${extra}/>`;
            };
            const label = (p, text, cls, dx = 0, dy = 0, anchor = "middle") => {
                const [x, y] = proj(p);
                return `<text x="${(x + dx).toFixed(1)}" y="${(y + dy).toFixed(1)}" text-anchor="${anchor}" class="${cls}">${text}</text>`;
            };
            // circle-arc in the plane spanned by orthonormal u, w, from u towards w
            const arcPts = (u, w, rad, ang, n = 48) => {
                const pts = [];
                for (let i = 0; i <= n; i++) {
                    const t = (ang * i) / n;
                    pts.push(P.add(P.scale(u, rad * Math.cos(t)), P.scale(w, rad * Math.sin(t))));
                }
                return pts;
            };
            const disc = (u, w, rad) => arcPts(u, w, rad, 2 * Math.PI, 96);

            const X = [1, 0, 0], Y = [0, 1, 0], Z = [0, 0, 1];
            const hHat = post.hmag > 1e-9 ? P.unit(post.h) : Z;
            const eHat = post.ecc > 1e-9 ? P.unit(post.evec) : null;
            let svg = "";

            // orbital planes as translucent discs
            svg += path(disc(X, Y, 1.75), "fill: var(--kg-plane-pre); stroke: var(--kg-rule)", 'stroke-width="1"');
            if (post.hmag > 1e-9) {
                const u = eHat || P.unit(cross0(hHat));
                svg += path(disc(u, P.cross(hHat, u), 1.75), "fill: var(--kg-plane-post); stroke: var(--kg-post); stroke-opacity: 0.35", 'stroke-width="1"');
            }

            // orbits
            const preOrbit = [];
            for (let i = 0; i <= 240; i++) {
                const E = (i / 240) * 2 * Math.PI;
                preOrbit.push([Math.cos(E) - pre.e, pre.q * Math.sin(E), 0]);
            }
            svg += path(preOrbit, "fill: none; stroke: var(--kg-ink)", 'stroke-width="1.5"');
            svg += path(conicPoints(post, 360, 3.5), "fill: none; stroke: var(--kg-post)",
                        'stroke-width="1.5" stroke-dasharray="6 4" clip-path="url(#kg-clip-3d)"');

            // original axes
            const AX = 2.05;
            svg += line([0, 0, 0], P.scale(X, AX), "stroke: var(--kg-muted)", 'stroke-width="1.1" marker-end="url(#kg-ah-ink)"');
            svg += line([0, 0, 0], P.scale(Y, AX), "stroke: var(--kg-muted)", 'stroke-width="1.1" marker-end="url(#kg-ah-ink)"');
            svg += line([0, 0, 0], P.scale(Z, 1.6), "stroke: var(--kg-muted)", 'stroke-width="1.1" marker-end="url(#kg-ah-ink)"');
            svg += label(P.scale(X, AX + 0.15), "x", "small");
            svg += label(P.scale(Y, AX + 0.15), "y", "small");
            // put the z and z' labels on opposite sides, on the side each arrow leans
            const zLeft = post.hmag > 1e-9 && proj(P.scale(hHat, 1.6))[0] >= proj(P.scale(Z, 1.6))[0];
            svg += label(P.scale(Z, 1.68), "z (h before)", "small", zLeft ? -6 : 6, 0, zLeft ? "end" : "start");

            // line of nodes (only defined when the plane tilts)
            let nodes = null;
            if (euler.kind === "general" || euler.kind === "disrupted") {
                const n = P.cross(Z, hHat);
                if (P.norm(n) > 1e-9) nodes = P.unit(n);
            }
            if (nodes) {
                svg += line(P.scale(nodes, -2.0), P.scale(nodes, 2.0), "stroke: var(--kg-ink)", 'stroke-width="1.2" stroke-dasharray="8 3 2 3"');
                // label whichever end of the line is lower in the view, away from z
                // (skipped when the line points at the viewer)
                const end = proj(nodes)[1] > proj(P.scale(nodes, -1))[1] ? 1 : -1;
                const [ex, ey] = proj(P.scale(nodes, 2.0 * end));
                if (Math.hypot(ex - CX, ey - CY) > 30) {
                    svg += label(P.scale(nodes, 2.0 * end), "line of nodes", "small", ex >= CX ? 4 : -4, 16,
                                 ex >= CX ? "start" : "end");
                }
            }

            // new axes: z' along h after, x' along the new periapsis
            if (post.hmag > 1e-9) {
                svg += line([0, 0, 0], P.scale(hHat, 1.6), "stroke: var(--kg-post)", 'stroke-width="2.2" marker-end="url(#kg-ah-post)"');
                svg += label(P.scale(hHat, 1.68), "z′ (h after)", "post", zLeft ? 6 : -6, 0, zLeft ? "start" : "end");
            }
            if (eHat && post.bound) {
                svg += line([0, 0, 0], P.scale(eHat, 1.45), "stroke: var(--kg-post)", 'stroke-width="2.2" marker-end="url(#kg-ah-post)"');
                svg += label(P.scale(eHat, 1.62), "x′ (new periapsis)", "post");
                if (showSecond.checked) {
                    svg += line([0, 0, 0], P.scale(eHat, -1.45), "stroke: var(--kg-post)", 'stroke-width="1.6" stroke-dasharray="2 3" marker-end="url(#kg-ah-post)"');
                    svg += label(P.scale(eHat, -1.62), "x-axis at 2nd SN (−x′)", "post");
                }
            }

            // Euler angle arcs
            const angleLabel = (u, w, rad, ang, text, cls) =>
                label(P.add(P.scale(u, rad * Math.cos(ang / 2)), P.scale(w, rad * Math.sin(ang / 2))), text, cls, 0, 4);
            const sub = '<tspan dy="3" font-size="9">E</tspan><tspan dy="-3"></tspan>';
            if (euler.thetaE > 1e-6 && euler.thetaE < Math.PI - 1e-6) {
                const w = P.unit(P.add(hHat, P.scale(Z, -hHat[2])));
                svg += path(arcPts(Z, w, 0.75, euler.thetaE), "fill: none; stroke: var(--kg-ink)", 'stroke-width="1.6"');
                svg += angleLabel(Z, w, 0.95, euler.thetaE, `θ${sub}`, "angle");
            }
            if (euler.kind === "general") {
                svg += path(arcPts(X, Y, 1.05, euler.phiE), "fill: none; stroke: var(--kg-vel)", 'stroke-width="1.8"');
                svg += angleLabel(X, Y, 1.25, euler.phiE, `φ${sub}`, "vel");
                if (euler.psiE !== null) {
                    const w = P.cross(hHat, nodes);
                    svg += path(arcPts(nodes, w, 0.85, euler.psiE), "fill: none; stroke: var(--kg-kick)", 'stroke-width="1.8"');
                    svg += angleLabel(nodes, w, 1.02, euler.psiE, `ψ${sub}`, "kick");
                }
            } else if ((euler.kind === "aligned" || euler.kind === "flipped") && eHat) {
                const ang = euler.kind === "aligned" ? euler.sumAngle : euler.diffAngle;
                const text = `φ${sub} ${euler.kind === "aligned" ? "+" : "−"} ψ${sub}`;
                svg += path(arcPts(X, Y, 1.05, ang), "fill: none; stroke: var(--kg-vel)", 'stroke-width="1.8"');
                svg += angleLabel(X, Y, 1.3, ang, text, "vel");
            }

            // the stars at the moment of the supernova
            const [qx, qy] = proj([0, 0, 0]);
            svg += `<circle cx="${qx.toFixed(1)}" cy="${qy.toFixed(1)}" r="6" style="fill: var(--kg-ink)"/>`;
            const [rx, ry] = proj(pre.r);
            svg += `<circle cx="${rx.toFixed(1)}" cy="${ry.toFixed(1)}" r="5" style="fill: var(--kg-panel); stroke: var(--kg-ink)" stroke-width="2"/>`;

            layer3d.innerHTML = svg;
        }

        // ---- readouts --------------------------------------------------------

        function renderReadouts(c) {
            const { s, pre, post, euler, vel } = c;
            set("ecc", s.e.toFixed(2));
            set("M", `${(+el.M.value).toFixed(0)}°`);
            set("theta", `${(+el.theta.value).toFixed(1)}°`);
            set("phi", `${(+el.phi.value).toFixed(0)}°`);

            set("loga", `${s.a < 10 ? s.a.toFixed(2) : s.a.toFixed(0)} R☉`);
            set("vk", `${s.vkkm.toFixed(0)} km/s`);
            set("m1", `${s.m1.toFixed(1)} M☉`);
            set("m2", `${s.m2.toFixed(1)} M☉`);
            set("dm", `${s.dm.toFixed(1)} M☉ (${(100 * s.dm / s.Mi).toFixed(0)}%)`);

            set("r-vc", `${s.vc.toFixed(0)} km/s`);
            // companion's speed about the pre-SN centre of mass (same frame as
            // the runaway speeds below): m1 / Mi of the relative orbital speed
            set("r-v2", `${(s.m1 / s.Mi * Math.hypot(pre.v[0], pre.v[1]) * s.vc).toFixed(0)} km/s`);
            set("r-t", `${(s.M / (2 * Math.PI)).toFixed(3)} P`);
            set("r-E", fmtAng(wrap360(pre.E / D)));
            set("r-nu", fmtAng(wrap360(pre.nu / D)));
            set("r-thM", fmtAng(wrap180((s.th - s.M) / D)));
            set("r-thnu", fmtAng(wrap180((s.th - pre.nu) / D)));
            set("r-thv", fmtAng(wrap180((s.th - pre.velAngle) / D)));
            set("r-dmcrit", `${(100 * pre.rmag / 2).toFixed(0)}% (r/2a)`);
            const st = out("r-state");
            // kick.f only treats e > 1 as disrupted; e = 1 happens with no kick
            // when exactly r/2a of the mass is lost
            const marginal = Math.abs(post.ecc - 1) < 1e-9;
            st.textContent = marginal ? "marginal (e = 1)" : post.bound ? "bound" : "unbound";
            st.className = "kg-state" + (post.bound && !marginal ? "" : " unbound");
            set("r-a", marginal ? "∞" : post.bound ? post.a.toFixed(3) : "—");
            set("r-e", post.ecc.toFixed(3));
            set("r-h", `${(post.hmag / pre.q).toFixed(3)}${post.h[2] < 0 ? ", reversed" : ""}`);
            const kms = v => `${(P.norm(v) * s.vc).toFixed(0)} km/s`;
            set("r-vcm", kms(vel.vcm));
            set("r-vsn", post.bound ? "—" : kms(vel.vExploding));
            set("r-vcomp", post.bound ? "—" : kms(vel.vCompanion));

            // Euler angles, as kick.f stores them in kick_info
            set("r-thetaE", fmtAng(euler.thetaE / D));
            set("r-phiE", euler.phiE === null ? "random" : fmtAng(euler.phiE / D));
            set("r-psiE", euler.psiE === null ? "random" : fmtAng(euler.psiE / D));
            const [tE, pE, sE] = ["θ<sub>E</sub>", "φ<sub>E</sub>", "ψ<sub>E</sub>"];
            const toPeri = a => (a === null ? "the angle to the new periapsis" : fmtAng(a / D));
            const notes = {
                general: `The orbit tilts by ${tE} about the line of nodes, which is ${pE} round from the x-axis. The new periapsis is ${sE} further round, in the new orbital plane.`,
                aligned: `The orbital plane does not tilt, so the line of nodes is undefined. COSMIC draws ${pE} at random and sets ${sE} so that ${pE} + ${sE} = ${toPeri(euler.sumAngle)}, the angle from the x-axis to the new periapsis.`,
                flipped: `The orbit reverses direction in the same plane, so the line of nodes is undefined. COSMIC draws ${pE} at random and sets ${sE} so that ${pE} − ${sE} = ${toPeri(euler.diffAngle)}, the angle from the x-axis to the new periapsis.`,
                disrupted: `The binary is disrupted, so only the tilt ${tE} of the (unbound) orbit is meaningful. COSMIC draws ${pE} and ${sE} at random.`,
            };
            let note = notes[euler.kind];
            if (euler.kind === "general" && euler.psiE === null) note += ` The new orbit is exactly circular, so ${sE} is random.`;
            out("r-euler-note").innerHTML = note;
        }

        function render() {
            const c = compute();
            renderTop(c);
            render3D(c);
            renderReadouts(c);
        }

        // ---- interaction -----------------------------------------------------

        // the remnant keeps at least 0.1 Msun, so the ejecta can't exceed m1 - 0.1
        function limitEjecta() {
            el.dm.max = (Math.round(10 * (+el.m1.value - 0.1)) / 10).toString();
            if (+el.dm.value > +el.dm.max) el.dm.value = el.dm.max;
        }
        el.m1.addEventListener("input", limitEjecta);
        limitEjecta();
        SLIDERS.forEach(id => el[id].addEventListener("input", render));
        showSecond.addEventListener("change", render);

        root.querySelectorAll("[data-kg-set]").forEach(b => b.addEventListener("click", () => {
            const s = state(), pre = P.preSupernova(s.e, s.M);
            const target = { M: s.M, nu: pre.nu, vel: pre.velAngle, in: pre.nu + Math.PI }[b.dataset.kgSet];
            el.theta.value = (Math.round(wrap360(target / D) * 10) / 10) % 360;
            render();
        }));

        root.querySelectorAll("[data-kg-view]").forEach(b => b.addEventListener("click", () => {
            const kind = b.dataset.kgView;
            if (kind === "reset") Object.assign(view, DEFAULT_VIEW);
            if (kind === "top") Object.assign(view, { az: -90 * D, el: 89.9 * D });
            if (kind === "nodes") {
                // look along the line of nodes, so the tilt appears edge-on
                const c = compute();
                const n = P.cross([0, 0, 1], c.post.h);
                const az = P.norm(n) > 1e-9 ? Math.atan2(n[1], n[0]) : 0;
                Object.assign(view, { az, el: 0 });
            }
            render();
        }));

        // drag (mouse, pen or touch) or arrow keys to rotate the 3D view
        let drag = null;
        svg3d.addEventListener("pointerdown", e => {
            drag = { x: e.clientX, y: e.clientY, az: view.az, el: view.el };
            svg3d.setPointerCapture(e.pointerId);
        });
        svg3d.addEventListener("pointermove", e => {
            if (!drag) return;
            view.az = drag.az - (e.clientX - drag.x) * 0.6 * D;
            view.el = Math.max(-89.9 * D, Math.min(89.9 * D, drag.el + (e.clientY - drag.y) * 0.6 * D));
            render();
        });
        const endDrag = () => { drag = null; };
        svg3d.addEventListener("pointerup", endDrag);
        svg3d.addEventListener("pointercancel", endDrag);
        svg3d.addEventListener("keydown", e => {
            const step = 5 * D;
            if (e.key === "ArrowLeft") view.az += step;
            else if (e.key === "ArrowRight") view.az -= step;
            else if (e.key === "ArrowUp") view.el = Math.min(89.9 * D, view.el + step);
            else if (e.key === "ArrowDown") view.el = Math.max(-89.9 * D, view.el - step);
            else return;
            e.preventDefault();
            render();
        });

        render();
    }

    function init() {
        document.querySelectorAll(".kick-geometry [data-kg-explorer]").forEach(initExplorer);
    }
    if (document.readyState === "loading") {
        document.addEventListener("DOMContentLoaded", init);
    } else {
        init();
    }
})();

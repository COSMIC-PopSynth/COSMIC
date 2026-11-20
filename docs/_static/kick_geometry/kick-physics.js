/*
    kick-physics.js

    The geometry of a supernova natal kick, following COSMIC's kick_pfahl
    (src/cosmic/src/kick.f) for a single supernova. This file has no DOM code
    so that it can be tested against the Fortran with node.

    Units: the pre-SN relative orbit has a = 1 and G M_i = 1, so speeds are in
    units of the circular speed V_c = (G M_i / a)^(1/2). After the supernova
    the total mass is mu = M_f / M_i (in units of M_i).

    Frame (as in kick.f): the pre-SN orbit lies in the x-y plane, +z is along
    the orbital angular momentum and +x points from the companion towards the
    periapsis of the exploding star's orbit about it.
*/
(function (root, factory) {
    const api = factory();
    if (typeof module === "object" && module.exports) {
        module.exports = api;
    } else {
        root.CosmicKickPhysics = api;
    }
})(typeof self !== "undefined" ? self : this, function () {
    "use strict";

    const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
    const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    const norm = a => Math.hypot(a[0], a[1], a[2]);
    const scale = (a, s) => a.map(v => v * s);
    const add = (a, b) => a.map((v, i) => v + b[i]);
    const unit = a => scale(a, 1 / norm(a));

    // angle between two vectors in [0, pi], as AngleBetweenVectors in kick.f
    function angleBetween(a, b) {
        return Math.acos(Math.max(-1, Math.min(1, dot(a, b) / (norm(a) * norm(b)))));
    }

    // Eccentric anomaly for mean anomaly M in [0, 2 pi] and 0 <= e < 1.
    // Newton-Raphson kept inside a bracket on the root (as in SolveKepler).
    function solveKepler(M, e) {
        const reflect = M > Math.PI;
        const Mr = reflect ? 2 * Math.PI - M : M;
        let lo = Mr, hi = Math.min(Mr + e, Math.PI);
        let E = Math.min(Math.max(Mr + 0.85 * e * (Mr > 0 ? 1 : 0), lo), hi);
        for (let i = 0; i < 100; i++) {
            const f = E - e * Math.sin(E) - Mr;
            if (Math.abs(f) < 1e-15) break;
            if (f < 0) lo = E; else hi = E;
            let En = E - f / (1 - e * Math.cos(E));
            if (!(En > lo && En < hi)) En = 0.5 * (lo + hi);
            if (En === E) break;
            E = En;
        }
        return reflect ? 2 * Math.PI - E : E;
    }

    // Position and velocity of the exploding star relative to its companion
    function preSupernova(e, M) {
        const E = solveKepler(M, e);
        const q = Math.sqrt(1 - e * e);
        const rmag = 1 - e * Math.cos(E);
        const r = [Math.cos(E) - e, q * Math.sin(E), 0];
        const v = [-Math.sin(E) / rmag, q * Math.cos(E) / rmag, 0];
        return {
            e, M, E, q, rmag, r, v,
            h: cross(r, v),
            nu: Math.atan2(r[1], r[0]),       // true anomaly
            velAngle: Math.atan2(v[1], v[0]), // direction of motion
        };
    }

    // natal kick with speed vk, azimuth theta and latitude phi (radians)
    function kickVector(vk, theta, phi) {
        return [vk * Math.cos(phi) * Math.cos(theta), vk * Math.cos(phi) * Math.sin(theta), vk * Math.sin(phi)];
    }

    // Orbit after the supernova: the separation is unchanged, the relative
    // velocity gains the kick and the total mass becomes mu (= M_f / M_i)
    function postSupernova(pre, kick, mu) {
        const v = add(pre.v, kick);
        const h = cross(pre.r, v);
        const hmag = norm(h);
        const evec = cross(v, h).map((c, i) => c / mu - pre.r[i] / pre.rmag);
        const ecc = norm(evec);
        const energy = dot(v, v) / 2 - mu / pre.rmag;
        return {
            v, h, hmag, evec, ecc, mu, energy,
            // kick.f treats the system as disrupted when e > 1
            bound: ecc <= 1,
            a: energy < 0 ? -mu / (2 * energy) : Infinity,
        };
    }

    /*
        Velocities of the stars after the supernova, in the pre-SN frame and
        relative to the pre-SN centre of mass, for an exploding star of mass
        m1 that leaves a remnant m1n and a companion m2 (as in kick.f). If the
        binary stays bound, vcm is its new systemic velocity. If it is
        disrupted, vExploding and vCompanion are the stars' velocities once
        they have escaped each other.
    */
    function velocities(pre, post, kick, m1, m2, m1n) {
        const Mi = m1 + m2, Mf = m1n + m2;
        const vcm = add(scale(pre.v, -m2 * (m1 - m1n) / (Mi * Mf)), scale(kick, m1n / Mf));
        if (post.bound) return { vcm };
        // relative velocity at infinity, from the eccentricity and angular
        // momentum vectors of the hyperbolic orbit
        const eHat = unit(post.evec), hHat = unit(post.h);
        const vinf = post.mu / post.hmag * Math.sqrt(post.ecc * post.ecc - 1);
        const vinfVec = add(scale(eHat, -vinf / post.ecc),
                            scale(cross(hHat, eHat), vinf * Math.sqrt(1 - 1 / (post.ecc * post.ecc))));
        return {
            vcm,
            vExploding: add(scale(vinfVec, m2 / Mf), vcm),
            vCompanion: add(scale(vinfVec, -m1n / Mf), vcm),
        };
    }

    /*
        Euler angles (z-x-z) of the post-SN orbital frame relative to the
        pre-SN frame, computed as in kick.f. A vector v' in the post-SN frame
        (x' = new periapsis, z' = new angular momentum) is
            v = Rz(phiE) Rx(thetaE) Rz(psiE) v'
        in the pre-SN frame. Angles that kick.f draws at random are returned
        as null, with `kind` saying why:
            "general"   - all three angles are defined
            "aligned"   - thetaE = 0, only phiE + psiE (= sumAngle) matters
            "flipped"   - thetaE = pi, only phiE - psiE (= diffAngle) matters
            "disrupted" - only thetaE is meaningful
        psiE is also null if the new orbit is exactly circular.
    */
    function eulerAngles(pre, post) {
        const thetaE = angleBetween(post.h, pre.h);
        const periAngle = post.ecc > 0 ? Math.atan2(post.evec[1], post.evec[0]) : null;
        if (!post.bound) {
            return { kind: "disrupted", thetaE, phiE: null, psiE: null };
        }
        if (thetaE === 0) {
            return { kind: "aligned", thetaE, phiE: null, psiE: null, sumAngle: periAngle };
        }
        if (thetaE === Math.PI) {
            return { kind: "flipped", thetaE, phiE: null, psiE: null, diffAngle: periAngle };
        }
        // line of nodes N = z x h, which is parallel to (-h_y, h_x, 0)
        const phiE = Math.atan2(post.h[0], -post.h[1]);
        let psiE = null;
        if (post.ecc > 0) {
            const nodes = cross(pre.h, post.h);
            const sign = dot(post.evec, pre.h) >= 0 ? 1 : -1;
            psiE = sign * angleBetween(post.evec, nodes);
        }
        return { kind: "general", thetaE, phiE, psiE };
    }

    // z-x-z rotation matrix R = Rz(phi) Rx(theta) Rz(psi), as ChangeBasis
    function rotationMatrix(thetaE, phiE, psiE) {
        const [ct, st, cp, sp, cs, ss] = [Math.cos(thetaE), Math.sin(thetaE), Math.cos(phiE),
                                          Math.sin(phiE), Math.cos(psiE), Math.sin(psiE)];
        return [
            [cp * cs - sp * ct * ss, -cp * ss - sp * ct * cs, st * sp],
            [sp * cs + cp * ct * ss, -sp * ss + cp * ct * cs, -st * cp],
            [st * ss, st * cs, ct],
        ];
    }

    return {
        cross, dot, norm, scale, add, unit, angleBetween,
        solveKepler, preSupernova, kickVector, postSupernova, velocities, eulerAngles, rotationMatrix,
    };
});

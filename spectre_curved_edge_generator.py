#!/usr/bin/python3

import numpy as np
import math
from time import time
import matplotlib.pyplot as plt

# some ideas from here: https://github.com/shrx/spectre
# and in particular for Python content: https://github.com/shrx/spectre/issues/1


class Spectre:

    # general-purpose linear interpolation function
    # super useful
    # x = argument
    # xa: input lower bound
    # xb: input upper bound
    # ya: output lower bound
    # yb: output upper bound

    def ntrp(x, xa, xb, ya, yb):
        return (x - xa) * (yb - ya) / (xb - xa) + ya

    PI = np.pi

    IDENTITY = np.array([[1, 0, 0], [0, 1, 0]], "float32")

    TILE_NAMES = [
        "Gamma",
        "Delta",
        "Theta",
        "Lambda",
        "Xi",
        "Pi",
        "Sigma",
        "Phi",
        "Psi",
    ]

    COLOR_MAP = {
        "Gamma": np.array((211, 95, 95), "f") / 255.0,
        "Gamma1": np.array((200, 55, 55), "f") / 255.0,
        "Gamma2": np.array((222, 135, 135), "f") / 255.0,
        "Delta": np.array((255, 179, 128), "f") / 255.0,
        "Theta": np.array((0, 255, 255), "f") / 255.0,
        "Lambda": np.array((128, 128, 128), "f") / 255.0,
        "Xi": np.array((145, 95, 211), "f") / 255.0,
        "Pi": np.array((95, 95, 211), "f") / 255.0,
        "Sigma": np.array((85, 212, 0), "f") / 255.0,
        "Phi": np.array((255, 170, 238), "f") / 255.0,
        "Psi": np.array((255, 221, 85), "f") / 255.0,
    }

    # generate the canonical spectre tile

    def generate_spectre_vertices():
        turtle_moves = (
            90,
            -60,
            90,
            60,
            -90,
            60,
            90,
            60,
            -90,
            60,
            0,
            60,
            90,
        )
        p = complex(0, 0)
        array = [(0.0, 0.0)]
        a = 90
        for m in turtle_moves:
            a += m
            ar = math.radians(a)
            p += complex(math.sin(ar), -math.cos(ar))
            array.append((p.imag, -p.real))
        return np.array(array, "float32")

    SPECTRE_POINTS = generate_spectre_vertices()

    SPECTRE_QUAD = SPECTRE_POINTS[[3, 5, 7, 11], :]

    def mul(A, B):
        AB = A.copy()
        AB[:, :2] = A[:, :2].dot(B[:, :2])
        AB[:, 2] += A[:, :2].dot(B[:, 2])
        return AB

    def __init__(
        self, iterations=1, amplitude=0.0, steps=32, image_w=19.2, image_h=10.8
    ):

        self.iterations = iterations
        self.amplitude = amplitude
        self.steps = steps
        self.image_w = image_w
        self.image_h = image_h

    # build_default_spectre()

    def buildSpectreBase(self):
        ttrans = np.array(
            [[1, 0, Spectre.SPECTRE_POINTS[8, 0]], [0, 1, Spectre.SPECTRE_POINTS[8, 1]]]
        )
        trot = np.array(
            [
                [np.cos(Spectre.PI / 6), -np.sin(Spectre.PI / 6), 0.0],
                [np.sin(Spectre.PI / 6), np.cos(Spectre.PI / 6), 0.0],
            ],
            "float32",
        )
        trsf = Spectre.mul(ttrans, trot)
        tiles = {}
        tiles["Gamma"] = MetaTile(
            tiles=[Tile("Gamma1"), Tile("Gamma2")],
            transformations=[Spectre.IDENTITY.copy(), trsf],
            quad=Spectre.SPECTRE_QUAD.copy(),
        )
        for label in Spectre.TILE_NAMES:
            if label != "Gamma":
                tiles[label] = Tile(label)
        return tiles

    def buildSupertiles(self, input_tiles):
        # First, use any of the nine-unit tiles in "tiles" to obtain a
        # list of transformation matrices for placing tiles within supertiles.
        quad = input_tiles["Delta"].quad

        transformations = [Spectre.IDENTITY.copy()]
        total_angle = 0
        trot = Spectre.IDENTITY.copy()
        transformed_quad = quad
        for _angle, _from, _to in (
            (Spectre.PI / 3, 3, 1),
            (0.0, 2, 0),
            (Spectre.PI / 3, 3, 1),
            (Spectre.PI / 3, 3, 1),
            (0.0, 2, 0),
            (Spectre.PI / 3, 3, 1),
            (-2 * Spectre.PI / 3, 3, 3),
        ):
            if _angle != 0:
                total_angle += _angle
                trot = np.array([[1, 0, 0], [0, 1, 0]]) * np.cos(
                    total_angle
                ) + np.array([[0, -1, 0], [1, 0, 0]]) * np.sin(total_angle)
                transformed_quad = quad.dot(trot[:, :2].T)  # + trot[:,2]
            last_trsf = transformations[-1]
            ttrans = Spectre.IDENTITY.copy()
            ttrans[:, 2] = (
                last_trsf[:, :2].dot(quad[_from, :])
                + last_trsf[:, 2]
                - transformed_quad[_to, :]
            )
            transformations.append(Spectre.mul(ttrans, trot))

        R = np.array([[-1, 0, 0], [0, 1, 0]], "float32")
        transformations = [Spectre.mul(R, trsf) for trsf in transformations]

        # Now build the actual supertiles, labeling appropriately.
        super_quad = quad[[2, 1, 2, 1], :]
        for i, itrsf in enumerate([6, 5, 3, 0]):
            trsf = transformations[itrsf]
            super_quad[i, :] = trsf[:, :2].dot(super_quad[i, :]) + trsf[:, 2]

        tiles = {}

        for label, substitutions in (
            ("Gamma", ("Pi", "Delta", None, "Theta", "Sigma", "Xi", "Phi", "Gamma")),
            ("Delta", ("Xi", "Delta", "Xi", "Phi", "Sigma", "Pi", "Phi", "Gamma")),
            ("Theta", ("Psi", "Delta", "Pi", "Phi", "Sigma", "Pi", "Phi", "Gamma")),
            ("Lambda", ("Psi", "Delta", "Xi", "Phi", "Sigma", "Pi", "Phi", "Gamma")),
            ("Xi", ("Psi", "Delta", "Pi", "Phi", "Sigma", "Psi", "Phi", "Gamma")),
            ("Pi", ("Psi", "Delta", "Xi", "Phi", "Sigma", "Psi", "Phi", "Gamma")),
            ("Sigma", ("Xi", "Delta", "Xi", "Phi", "Sigma", "Pi", "Lambda", "Gamma")),
            ("Phi", ("Psi", "Delta", "Psi", "Phi", "Sigma", "Pi", "Phi", "Gamma")),
            ("Psi", ("Psi", "Delta", "Psi", "Phi", "Sigma", "Psi", "Phi", "Gamma")),
        ):
            tiles[label] = MetaTile(
                tiles=[input_tiles[subst] for subst in substitutions if subst],
                transformations=[
                    trsf for subst, trsf in zip(substitutions, transformations) if subst
                ],
                quad=super_quad,
            )
        return tiles

    def render(self):
        start = time()
        tiles = self.buildSpectreBase()
        for _ in range(self.iterations):
            tiles = self.buildSupertiles(tiles)
        time1 = time() - start
        print(f"supertiling loop: {time1:.04f} seconds")

        start = time()
        polygons = []
        tiles["Delta"].draw(polygons)
        plt.style.use("dark_background")
        if self.image_w > 0 and self.image_h > 0:
            plt.rcParams["figure.figsize"] = (self.image_w, self.image_h)
        plt.axis("equal")
        plt.tight_layout(pad=1)
        for pts, color in polygons:
            self.apply_edge_sine_and_plot(pts[:, 1], pts[:, 0], color)
        time2 = time() - start
        print(
            f"tile recursion loop and image generation: {time2:.04f} seconds, generated {len(polygons)} tiles"
        )
        if self.image_w > 0 and self.image_h > 0:
            plt.savefig(f"output_image_{self.image_w}w_{self.image_h}h.png")
        plt.show()

    def apply_edge_sine_and_plot(self, nxp, nyp, color):
        xp = nxp.tolist()
        yp = nyp.tolist()
        # to close rendered polygons,
        # must wrap these arrays around
        xp += [xp[0]]
        yp += [yp[0]]
        # wiggle spectre tile edges if specified
        if self.steps > 1 and self.amplitude > 0:
            xs = []
            ys = []
            # internal normalization of user-entered amplitude
            amp = self.amplitude * 0.3333
            for p in range(0, len(xp) - 1):
                # xa,xb and ya,yb are the cartesian
                # vertex coordinates for a spectre edge
                xa, xb = xp[p], xp[p + 1]
                ya, yb = yp[p], yp[p + 1]
                # convert edge cartesian to polar coordinates
                angle = math.atan2((yb - ya), (xb - xa))
                radius = math.sqrt((yb - ya) ** 2 + (xb - xa) ** 2)
                for n in range(self.steps):
                    # create wiggle waveform using sine function
                    arg = Spectre.ntrp(n, 0, self.steps, 0, math.pi * 2)
                    sw = -math.sin(arg) * radius * amp
                    # mix orthogonal edge and sinewave components
                    rr = Spectre.ntrp(n, 0, self.steps, 0, radius)
                    xx = rr * math.cos(angle) + sw * math.sin(angle)
                    yy = rr * math.sin(angle) - sw * math.cos(angle)
                    xs += [xx + xa]
                    ys += [yy + ya]
            xp = xs
            yp = ys

        plt.fill(xp, yp, facecolor=color)  # color)
        # plot dark border around each tile
        plt.plot(xp, yp, color="black", linewidth=0.25)
        # some debugging code
        # for n in range(len(xp)):
        #    plt.plot(xp[n], yp[n], marker="o", label=f"{n}")
        #    plt.text(xp[n], yp[n], f"{n}")
        # plt.legend()


class Tile:
    def __init__(self, label):
        self.label = label
        self.quad = Spectre.SPECTRE_QUAD.copy()

    def draw(self, polygons, tile_transformation=Spectre.IDENTITY.copy()):
        vertices = (
            # build_default_spectre()
            Spectre.SPECTRE_POINTS.dot(tile_transformation[:, :2].T)
            + tile_transformation[:, 2]
        )
        # print(vertices)
        # quit()
        polygons.append((vertices, Spectre.COLOR_MAP[self.label]))


class MetaTile:
    def __init__(self, tiles=[], transformations=[], quad=Spectre.SPECTRE_QUAD.copy()):
        self.tiles = tiles
        self.transformations = transformations
        self.quad = quad

    def draw(self, polygons, transformation=Spectre.IDENTITY.copy()):
        for tile, trsf in zip(self.tiles, self.transformations):
            tile.draw(polygons, Spectre.mul(transformation, trsf))


# class arguments:
# iterations controls recursions
# iterations determines recursive depth:
# iterations tiles
# 0            1
# 1            9
# 2           71
# 3          559
# 4         4401
# 5         Don't even think about it
# steps determines the detail of the edge waveforms
# amplitude controls the depth of the edge waveform curvature
# 0 = no curvature
# 1 = maximum curvature consistent with interlocking tiles
# image_w, image_h determine the size of an exported image
# using matplotlib's image size convention


def main():
    spectre = Spectre(iterations=1, steps=32, amplitude=1, image_w=19.2, image_h=10.8)

    spectre.render()


if __name__ == "__main__":
    main()

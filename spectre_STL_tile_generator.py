#!/usr/bin/env python
# type: ignore

import matplotlib.pyplot as plt
import cadquery as cq
import math, os, sys

# class arguments:

# steps: how many discrete sections each edge wiggle has

# tiles: how many graphic images and STL files to create

# edge_length (mm): controls the overall
# size of a tile, which is derived from this length

# amplitude: depth of edge waveform, an excessive
# value may prevent tiles from linking together

# shrink_factor: shrinks the tile size
# center-to-edge, units mm,
# to allow interlocking of printed tiles 
# without too much force, an excessive
# value may prevent creation of printable STL files

# data_dir is an automatically created
# local subdirectory to receive
# graphic images and STL files 

class PrintSpectre:
    def __init__(
        self,
        steps = 32,
        tiles = 8,
        edge_length=12,
        amplitude=1.0,
        shrink_factor=0.25,
        data_dir="spectre_tiles",
    ):
        self.steps = steps
        self.tiles = tiles
        self.edge_length = edge_length
        self.amplitude = amplitude
        self.shrink_factor = shrink_factor
        self.data_dir = data_dir

    # linear interpolation function
    def ntrp(self, x, xa, xb, ya, yb):
        return (x - xa) * (yb - ya) / (xb - xa) + ya
    
    # generate the canonical spectre tile

    def gen_spectre(self):
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
        xp, yp = [0], [0]
        a = 180
        for m in turtle_moves:
            a += m
            ar = math.radians(a)
            p += self.edge_length * complex(math.sin(ar), -math.cos(ar))
            xp += [-p.imag]
            yp += [p.real]
        xp += [0]
        yp += [0]
        return xp, yp

    # a two-pass perimeter shrinking scheme
    # to allow parts to fit snugly without jamming
    # this function is called twice:
    # cw, then ccw
    def shrink_perimeter_core(self, xp, yp, sign=1):
        if sign > 0:
            xp = xp[::-1]
            yp = yp[::-1]
        xs = []
        ys = []
        for n in range(len(xp) - 1):
            xa = xp[n]
            xb = xp[n + 1]
            ya = yp[n]
            yb = yp[n + 1]
            dx = xb - xa
            dy = yb - ya
            angle = math.atan2(dy, dx)
            sx = -math.sin(angle) * self.shrink_factor * sign * 0.5
            sy = math.cos(angle) * self.shrink_factor * sign * 0.5
            xs += [-sx + xa]
            ys += [-sy + ya]
        if sign > 0:
            xs = xs[::-1]
            ys = ys[::-1]
        return xs, ys

    # we call the perimeter shrink function
    # twice, once CW aound the tile, once CCW,
    #  then average the results
    def shrink_perimeter(self, xp, yp):
        # resx = []
        # resy = []
        # Each call to the shrink_perimeter_core
        # function trims the specified factor
        # so apply 1/2 the specified factor, twice
        # traverse clockwise
        xs1, ys1 = self.shrink_perimeter_core(xp, yp, -1)
        # traverse counterclockwise
        xs2, ys2 = self.shrink_perimeter_core(xp,yp, 1)
        # combine results
        resx = [(a + b) * .5 for a, b in zip(xs1, xs2)]
        resy = [(a + b) * .5 for a, b in zip(ys1, ys2)]
        return resx,resy
        
    # xs,ys, cartesian Spectre tile coordinates
    # amplitude 0 = classic spectre with straight edges
    # aplitude 1.0 = extreme spectre with deep wiggles

    def build_spectre(self, xs, ys, amplitude, plot=False):
        print(f"Building spectre with edge amplitude {amplitude:.02f} ...")
        sys.stdout.flush()
        # internal amplitude conversion
        amp = amplitude * 0.3333
        xp = []
        yp = []
        for p in range(len(xs) - 1):
            # xa,xb and ya,yb are the
            # cartesian vertices for a spectre edge
            xa, xb = xs[p], xs[p + 1]
            ya, yb = ys[p], ys[p + 1]
            # convert edge to polar coordinates
            angle = math.atan2((yb - ya), (xb - xa))
            radius = math.sqrt((yb - ya) ** 2 + (xb - xa) ** 2)
            for n in range(self.steps):
                # create wiggle waveform using sine function
                arg = self.ntrp(n, 0, self.steps, 0, math.pi * 2)
                sw = math.sin(arg) * radius * amp
                # mix orthogonal edge and sinewave components
                rr = self.ntrp(n, 0, self.steps, 0, radius)
                xx = rr * math.cos(angle) + sw * math.sin(angle)
                yy = rr * math.sin(angle) - sw * math.cos(angle)
                xp += [xx + xa]
                yp += [yy + ya]
        # add one more data point
        xp += [0]
        yp += [0]

        # shrink the object's perimeter to avoid
        # an overly tight fit when printed and assembled
        # this is essential for wiggly edges
        xs, ys = self.shrink_perimeter(xp, yp)

        # create STL files for each version
        # "thickness_mm" is the STL tile thickness in mm
        thickness_mm = 2
        points = [[x, y] for x, y in zip(xs, ys)]

        if True: # create STL files
            # create cadquery object for spectre tile
            spectre = cq.Workplane("XY").polyline(points).close().extrude(thickness_mm)

            # save the spectre object in STL format
            cq.exporters.export(
                spectre, f"{self.data_dir}/spectre_tile_amplitude_{amplitude:.01f}.stl"
            )

        if plot:
            # plot and save spectre tile images
            plt.close()
            plt.style.use("dark_background")
            fig, ax = plt.subplots(figsize=(19.2, 10.8))
            ax.set_aspect("equal")
            plt.xlabel("Size mm", fontsize=14)
            plt.ylabel("Size mm", fontsize=14)
            plt.fill(xp, yp, color="#004000")
            plt.plot(xp, yp, color="red", label="Generated size")
            plt.plot(
                xs,
                ys,
                color="#4080ff",
                label=f"Shrunk by {self.shrink_factor:0.3f} mm",
            )

            plt.legend(fontsize=16)
            plt.title(f"Spectre tile with edge amplitude {amplitude:.02f}", fontsize=20)
            plt.savefig(f"{self.data_dir}/spectre_image_amplitude_{amplitude:.02f}.png")
            # plt.show()


def main():
    data_dir = "spectre_tiles"
    if not os.path.isdir(data_dir):
        os.mkdir(data_dir)
    # steps in the edge wiggle waveform 
    steps = 32
    # number of versions to create
    tiles = 8
    edge_length = 12 # mm per spectre edge
    amplitude = 1 # depth of wiggles
    # amount of shrinkage to allow
    #  parts to easily fit together
    shrink_factor = 0.25 # units mm
    # default spectre edge length mm
    edge_length = 12
    spectre = PrintSpectre(steps,tiles,edge_length, amplitude, shrink_factor, data_dir)
    # create the basic spectre data set
    xs, ys = spectre.gen_spectre()
    spectre.spx = xs
    spectre.spy = ys
    # tiles = number of tile versions moving 
    # between flat sides and deep wiggles
    for n in range(spectre.tiles):
        # apply edge wiggle amplitudes between 0 and 1
        amp = spectre.ntrp(n, 0, tiles - 1, 0, 1)
        spectre.build_spectre(xs, ys, amp, True)
    print(f"Done, created {tiles} Spectre tiles and images.")


if __name__ == "__main__":
    main()

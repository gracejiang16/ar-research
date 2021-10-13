# Python file to calculate points for the twoxsix lower scheme.
# Author: Grace Jiang
# Date: July 2021

import math

import numpy
import numpy as np

import generic
import matplotlib

from numpy import matrix
from numpy import linalg


# data1a is desired resulting name of file for first set of data for the "one" annulus, data1b is for "three" annulus
# this method exhaustively searches for all points on the circle around each antipode
# two "chunks" of code here because this scheme isn't symmetrical (2x6)
def exhaust_two(data1_one):
    open(data1_one, "w").close()  # clears file if it exists todo

    distance_small = 73.5 / 2  # degrees
    distance_big = 94.5 / 2

    big_radius = (360 - 4 * distance_small) / 2
    small_radius = (360 - 4 * distance_big) / 2

    grid_increment = 0.8  # degrees
    max_error = 0.4  # degrees

    far_lat_min = -60
    far_lat_max = 60
    near_lat_min = -60
    near_lat_max = 60

    far_lon_min = 0
    far_lon_max = 0.1  # don't need lon loop
    near_lon_min = 90
    near_lon_max = 270

    far_current_lat = far_lat_min  # starts at min, goes to max
    far_current_lon = far_lon_min

    f = open(data1_one, "w")  # data1a (and b) comes from this exhaustive search; lat changes but lon is constant

    while far_current_lat <= far_lat_max:
        far_current_lon = far_lon_min

        lat_2 = near_lat_min  # degrees (these are each grid square on the near side)
        lon_2 = near_lon_min  # degrees

        while (far_current_lon <= far_lon_max):
            lat_1 = far_current_lat * -1  # degrees (finding antipode)
            lon_1 = far_current_lon + 180  # degrees (finding antipode)

            # go through all different radii
            current_radius = small_radius
            while current_radius <= big_radius:
                lat_2 = near_lat_min
                while lat_2 <= near_lat_max:
                    lon_2 = near_lon_min
                    while lon_2 <= near_lon_max:

                        # delta is the calculated arc length between antipode and every near side point todo changed arclength
                        delta = math.acos((math.sin(math.radians(lat_1))) * (math.sin(math.radians(lat_2)))
                                          + (math.cos(math.radians(lat_1))) * (math.cos(math.radians(lat_2))) *
                                          (math.cos(math.radians(lon_1 - lon_2))))  # radians

                        if abs(delta - math.radians(current_radius)) <= math.radians(max_error):
                            s = str(lat_2) + " " + str(lon_2) + " "  # lat, lon
                            f.write(s + "\n")

                        lon_2 += grid_increment
                    lat_2 += grid_increment

                line_separator = "far lat: " + str(
                    far_current_lat) + " far lon: " + str(
                    far_current_lon) + " radius: " + str(
                    current_radius) + "\n"  # keep track of diff groups of points for diff far side points
                f.write(line_separator)

                current_radius += 0.8

            far_current_lon += 6  # lon interval = 6 todo

        far_current_lat += grid_increment

    f.close()
    print("--- end of exhaust: two ---")  # todo


def exhaust_six(data1_three):
    open(data1_three, "w").close()  # clears file if it exists todo
    # second chunk of code for the "three" annulus:

    distance_small = 220.5 / 6  # "one" annulus distance (degrees)
    distance_big = 283.5 / 6

    big_radius = (12 * distance_big - 360) / 2  # 103.5
    small_radius = (12 * distance_small - 360) / 2  # 40.5

    grid_increment = 0.8  # degrees
    max_error = 0.4  # degrees

    far_lat_min = -60
    far_lat_max = 60
    near_lat_min = -60
    near_lat_max = 60

    far_lon_min = 0
    far_lon_max = 0.1  # don't need lon loop
    near_lon_min = 90
    near_lon_max = 270

    far_current_lat = far_lat_min  # starts at min, goes to max
    far_current_lon = far_lon_min

    # second chunk for "three" annulus

    ff = open(data1_three, "w")

    while far_current_lat <= far_lat_max:
        far_current_lon = far_lon_min

        lat_2 = near_lat_min  # degrees (these are each grid square on the near side)
        lon_2 = near_lon_min  # degrees

        while (far_current_lon <= far_lon_max):
            lat_1 = far_current_lat * -1  # degrees (finding antipode)
            lon_1 = far_current_lon + 180  # degrees (finding antipode)

            # go through all different radii
            current_radius = small_radius
            while current_radius <= big_radius:
                lat_2 = near_lat_min
                while lat_2 <= near_lat_max:
                    lon_2 = near_lon_min
                    while lon_2 <= near_lon_max:

                        # delta is the calculated arc length between antipode and every near side point todo changed arc length
                        delta = math.acos((math.sin(math.radians(lat_1))) * (math.sin(math.radians(lat_2)))
                                          + (math.cos(math.radians(lat_1))) * (math.cos(math.radians(lat_2))) *
                                          (math.cos(math.radians(lon_1 - lon_2))))  # radians

                        if abs(delta - math.radians(current_radius)) <= math.radians(max_error):
                            midpoint = get_midpoint(far_current_lat, far_current_lon, lat_2, lon_2)
                            s = str(lat_2) + " " + str(lon_2) + " " + \
                                str(midpoint[0]) + " " + str(midpoint[1])  # lat, lon, midlat, midlon
                            ff.write(s + "\n")

                        lon_2 += grid_increment
                    lat_2 += grid_increment

                line_separator = "far lat: " + str(
                    far_current_lat) + " far lon: " + str(
                    far_current_lon) + " radius: " + str(
                    current_radius) + "\n"  # keep track of diff groups of points for diff far side points
                ff.write(line_separator)

                current_radius += 0.8

            far_current_lon += 18  # lon interval = 18 degrees; computational complexity is too much for todo
            # the 2x6 lower scheme; so must  select 18-degree as the interval todo

        far_current_lat += grid_increment

    ff.close()


def remove_repeats_two(coors_one, resulting_no_repeats_one, max_diff):
    max_difference = max_diff

    open(resulting_no_repeats_one, 'w').close()  # to clear all present text todo

    f = open(coors_one, 'r')

    lats = []
    lons = []  # originals (reset for each group of points)

    for line in f:
        line = line.strip()

        if "l" not in line:  # sectioning for diff groups of points
            c1, c2 = line.split(' ')
            lats.append(float(c1))
            lons.append(float(c2))

        elif "l" in line:  # sectioning for diff groups of points

            i = 0
            while i < len(lats):
                j = i + 1
                while j < len(lats):
                    if generic.distance(lats[i], lons[i], lats[j], lons[j]) <= math.radians(max_difference):
                        # print(lats[j], lons[j])
                        del lats[j]
                        del lons[j]
                    else:
                        j += 1
                i += 1

            ff = open(resulting_no_repeats_one, 'a')
            for lat, lon in zip(lats, lons):
                ff.write(str(lat) + " " + str(lon) + " " + "\n")

            ff.write(line)
            ff.write("\n")
            ff.close()

            lats = []
            lons = []
    f.close()
    print("--- end of remove_repeats_two ---")  # todo


# removes points from coors that are too close together according to the "max_diff" variable
def remove_repeats_six(coors_three, resulting_no_repeats_three, max_diff):
    max_difference = max_diff

    open(resulting_no_repeats_three, 'w').close()  # to clear all present text

    f = open(coors_three, 'r')

    lats = []
    lons = []  # originals (reset for each group of points)

    midlats = []
    midlons = []  # to retain midpoints

    for line in f:
        line = line.strip()

        if "l" not in line:  # sectioning for diff groups of points
            c1, c2, c3, c4 = line.split(' ')
            lats.append(float(c1))
            lons.append(float(c2))
            midlats.append(float(c3))
            midlons.append(float(c4))

        elif "l" in line:  # sectioning for diff groups of points

            i = 0
            while i < len(lats):
                j = i + 1
                while j < len(lats):
                    if generic.distance(lats[i], lons[i], lats[j], lons[j]) <= math.radians(max_difference):
                        # print(lats[j], lons[j])
                        del lats[j]
                        del lons[j]
                        del midlats[j]
                        del midlons[j]
                    else:
                        j += 1
                i += 1

            ff = open(resulting_no_repeats_three, 'a')
            for lat, lon, midlat, midlon in zip(lats, lons, midlats, midlons):
                ff.write(str(lat) + " " + str(lon) + " " + str(midlat) + " " + str(midlon) + "\n")

            ff.write(line)
            ff.write("\n")
            ff.close()

            lats = []
            lons = []
            midlats = []
            midlons = []  # to reset
    f.close()


# finds the middle skip coordinates from the 2nd skip "coors" and puts them in resulting_corresponding
def find_middle_skip(coors, resulting_corresponding):
    open(resulting_corresponding, 'w').close()  # to clear all present text

    f = open(coors, 'r')

    lats = []
    lons = []  # originals (reset for each group of points)

    new_lats = []
    # new_lons = []  # points we're looking for (reset for each group of points) todo comment this out

    for line in f:
        line = line.strip()

        if "l" not in line:  # sectioning for diff groups of points
            lat, lon = line.split(' ')
            lats.append(float(lat))
            lons.append(float(lon))

        elif "l" in line:  # sectioning for diff groups of points
            fsp = get_coors_from_str(line)  # extract far side point from line

            lat_fsp = fsp[0]
            lon_fsp = fsp[1]

            for lati, long in zip(lats, lons):  # actual finding corres point with helper method get_corres todo woth=with
                current_corres = get_mid_corres_normal(lat_fsp, lon_fsp, lati, long)
                corres_lat = current_corres
                new_lats.append(corres_lat)
                # corres_lon = current_corres[1]  # we don't need this part ONLY IF we disregard the lons from get_corres()
                # new_lons.append(corres_lon)

            ff = open(resulting_corresponding, 'a')  # writing corres points into resulting_corresponding
            # for lat, lon in zip(new_lats, new_lons):
            #     ff.write(str(lat) + " " + str(lon) + "\n")

            for lat in zip(new_lats):
                ff.write(
                    str(lat[0]) + "\n")  # "lat" is a tuple, so just access the first (and only) element
                # convert the lat into degrees, then to string, then write onto file

            ff.write(line)
            ff.write("\n")
            ff.close()

            lats = []
            lons = []  # to reset

            new_lats = []
            new_lons = []  # to reset

    f.close()
    print("--- end of find_middle_skip ---")  # todo


# finds the second skip from third skip file coors_three and puts results in resulting_corresponding_second_skip
def find_fourth_skip(coors_three, resulting_corresponding_second_skip):
    open(resulting_corresponding_second_skip, 'w').close()  # to clear all present text

    f = open(coors_three, 'r')

    lats = []
    lons = []
    midlats = []
    midlons = []  # originals (reset for each group of points)

    new_lats = []
    new_lons = []  # points we're looking for (reset for each group of points)

    for line in f:
        line = line.strip()

        if "l" not in line:  # sectioning for diff groups of points
            lat, lon, midlat, midlon = line.split(' ')
            lats.append(float(lat))
            lons.append(float(lon))
            midlats.append(float(midlat))
            midlons.append(float(midlon))

        elif "l" in line:  # sectioning for diff groups of points
            fsp = get_coors_from_str(line)  # extract far side point and alpha from line
            lat_fsp = fsp[0]
            lon_fsp = fsp[1]
            alpha = fsp[2] + 180  # need to +180 since read in "radius" only  # todo

            # print(lat_fsp, lon_fsp, alpha)
            for lati, long, midlati, midlong in zip(lats,
                                                    lons, midlats,
                                                    midlons):  # actual finding corres point (second skip) with helper method get_second_skip

                if midlong < -90:
                    midlong += 360

                second_skip = get_second_skip(lati, long, midlati, midlong,
                                              alpha / 2)  # finding antipodes of midpoints, new alpha is half its original
                second_lat = second_skip
                new_lats.append(float(second_lat))  # recording the second skip lat

                # second_lon = second_skip[1]
                # new_lons.append(second_lon)  # recording the second skip lon (no need)

            ff = open(resulting_corresponding_second_skip,
                      'a')  # writing corres points into resulting_corresponding_second_skip
            # for lat, lon in zip(new_lats, new_lons):
            #     ff.write(str(lat) + " " + str(lon) + "\n")

            for lat in zip(new_lats):
                # ff.write(str(math.degrees(lat)) + " " +
                #          str(math.degrees(lon)) + "\n")  # "lat" and "lon" are tuples, so access the first (and only) element of each
                # convert the lat and lon into degrees, then to string, then write onto file
                # print(lat, math.degrees(lat[0]))
                ff.write(str(math.degrees(lat[0])) + "\n")

            ff.write(line + "\n")
            ff.close()

            lats = []
            lons = []
            midlats = []
            midlons = []  # to reset

            new_lats = []
            new_lons = []  # to reset

    f.close()
    print("--- end of find_fourth_skip ---")  # todo


# finds the 5th skip from third skip file coors_three and puts results in resulting_corresponding_second_skip
def find_fifth_skip(coors_three, resulting_corresponding_second_skip):
    open(resulting_corresponding_second_skip, 'w').close()  # to clear all present text

    f = open(coors_three, 'r')

    lats = []
    lons = []
    midlats = []
    midlons = []  # originals (reset for each group of points)

    new_lats = []
    new_lons = []  # points we're looking for (reset for each group of points)

    for line in f:
        line = line.strip()

        if "l" not in line:  # sectioning for diff groups of points
            lat, lon, midlat, midlon = line.split(' ')
            lats.append(float(lat))
            lons.append(float(lon))
            midlats.append(float(midlat))
            midlons.append(float(midlon))

        elif "l" in line:  # sectioning for diff groups of points
            fsp = get_coors_from_str(line)  # extract far side point and alpha from line
            lat_fsp = fsp[0]
            lon_fsp = fsp[1]
            alpha = fsp[2] + 180  # need to +180 since read in "radius" only todo

            # print(lat_fsp, lon_fsp, alpha)
            for lati, long, midlati, midlong in zip(lats,
                                                    lons, midlats,
                                                    midlons):  # actual finding corres point (second skip) with helper method get_second_skip

                if midlong < -90:
                    midlong += 360

                second_skip = get_second_skip(midlati, midlong, lati, long,
                                              alpha / 2)  # finding antipodes of midpoints, new alpha is half its original todo where do we find the antipode??
                second_lat = second_skip
                new_lats.append(float(second_lat))  # recording the second skip lat

                # second_lon = second_skip[1]
                # new_lons.append(second_lon)  # recording the second skip lon (no need)

            ff = open(resulting_corresponding_second_skip, 'a')  # writing corres points into resulting_corresponding_second_skip
            # for lat, lon in zip(new_lats, new_lons):
            #     ff.write(str(lat) + " " + str(lon) + "\n")

            for lat in zip(new_lats):
                ff.write(str(math.degrees(lat[0])) + "\n")

            ff.write(line + "\n")
            ff.close()

            lats = []
            lons = []
            midlats = []
            midlons = []  # to reset

            new_lats = []
            new_lons = []  # to reset

    f.close()
    print("--- end of find-5th-skip ---")


def find_sixth_skip(coors_three, resulting_corresponding_second_skip):
    open(resulting_corresponding_second_skip, 'w').close()  # to clear all present text

    f = open(coors_three, 'r')

    lats = []
    lons = []
    midlats = []
    midlons = []  # originals (reset for each group of points)

    new_lats = []
    new_lons = []  # points we're looking for (reset for each group of points)

    for line in f:
        line = line.strip()

        if "l" not in line:  # sectioning for diff groups of points
            lat, lon, midlat, midlon = line.split(' ')
            lats.append(float(lat))
            lons.append(float(lon))
            midlats.append(float(midlat))
            midlons.append(float(midlon))

        elif "l" in line:  # sectioning for diff groups of points
            fsp = get_coors_from_str(line)  # extract far side point and alpha from line
            lat_fsp = fsp[0]
            lon_fsp = fsp[1]
            alpha = fsp[2] + 180  # Need to +180 since read in 'radius' only

            # print(lat_fsp, lon_fsp, alpha)
            for lati, long, midlati, midlong in zip(lats,
                                                    lons, midlats,
                                                    midlons):  # actual finding corres point (second skip) with helper method get_second_skip

                if midlong < -90:
                    midlong += 360

                # second_skip = get_second_skip(midlati, midlong, lati, long,
                #                               alpha / 2)  # finding antipodes of midpoints, new alpha is half its original
                # second_lat = second_skip
                second_lat = lati
                new_lats.append(float(second_lat))  # recording the second skip lat

                # second_lon = second_skip[1]
                # new_lons.append(second_lon)  # recording the second skip lon (no need)

            ff = open(resulting_corresponding_second_skip,
                      'a')  # writing corres points into resulting_corresponding_second_skip
            # for lat, lon in zip(new_lats, new_lons):
            #     ff.write(str(lat) + " " + str(lon) + "\n")

            for lat in zip(new_lats):
                # ff.write(str(math.degrees(lat)) + " " +
                #          str(math.degrees(lon)) + "\n")  # "lat" and "lon" are tuples, so access the first (and only) element of each
                # convert the lat and lon into degrees, then to string, then write onto file
                # print(lat, math.degrees(lat[0]))
                ff.write(str(lat[0]) + "\n")

            ff.write(line + "\n")
            ff.close()

            lats = []
            lons = []
            midlats = []
            midlons = []  # to reset

            new_lats = []
            new_lons = []  # to reset

    f.close()
    print("--- end of find_sixth_skip ---")


# finds first skip point in the "three" annulus by cutting half from the fsp and the newly found second skip point
def find_first_skip(coors_three, resulting_corresponding_first_skip):
    open(resulting_corresponding_first_skip, 'w').close()  # to clear all present text

    f = open(coors_three, 'r')

    lats = []
    lons = []
    midlats = []  # we use mid-antipodes as v2 and fsp as v1
    midlons = []  # originals (reset for each group of points)

    new_lats = []
    new_lons = []  # points we're looking for (reset for each group of points)

    for line in f:
        line = line.strip()

        if "l" not in line:  # sectioning for diff groups of points
            lat, lon, midlat, midlon = line.split(' ')
            lats.append(float(lat))
            lons.append(float(lon))
            midlats.append(float(midlat))
            midlons.append(float(midlon))

        elif "l" in line:  # sectioning for diff groups of points
            fsp = get_coors_from_str(line)  # extract far side point from line (and radius)
            lat_fsp = fsp[0]
            lon_fsp = fsp[1]
            alpha = fsp[2] + 180

            # print(lat_fsp,lon_fsp,alpha)
            for midlati, midlong in zip(midlats,
                                        midlons):  # actual finding corres point (first skip) with helper method get_first_skip
                first_skip = get_first_skip(midlati, midlong, lat_fsp, lon_fsp, alpha / 2)
                first_lat = first_skip
                new_lats.append(first_lat)
                # first_lon = first_skip[1]  # don't include lons
                # new_lons.append(first_lon)

            ff = open(resulting_corresponding_first_skip,
                      'a')  # writing corres points into resulting_corresponding_first_skip
            # for lat, lon in zip(new_lats, new_lons):
            #     ff.write(str(lat) + " " + str(lon) + "\n")

            for lat in zip(new_lats):
                # print(math.degrees(lat[0]))
                ff.write(str(math.degrees(
                    lat[0])) + "\n")  # "lat" and "lon" are tuples, so access the first (and only) element of each
                # convert the lat and lon into degrees, then to string, then write onto file

            ff.write(line)
            ff.write("\n")
            ff.close()

            lats = []
            lons = []  # to reset
            midlats = []
            midlons = []  # to reset

            new_lats = []
            new_lons = []  # to reset

    f.close()
    print("--- end of find_first_skip ---")


def find_second_skip(coors_three, resulting_corresponding_first_skip):
    open(resulting_corresponding_first_skip, 'w').close()  # to clear all present text

    f = open(coors_three, 'r')

    lats = []
    lons = []
    midlats = []  # we use mid-antipodes as v2 and fsp as v1
    midlons = []  # originals (reset for each group of points)

    new_lats = []
    new_lons = []  # points we're looking for (reset for each group of points)

    for line in f:
        line = line.strip()

        if "l" not in line:  # sectioning for diff groups of points
            lat, lon, midlat, midlon = line.split(' ')
            lats.append(float(lat))
            lons.append(float(lon))
            midlats.append(float(midlat))
            midlons.append(float(midlon))

        elif "l" in line:  # sectioning for diff groups of points
            fsp = get_coors_from_str(line)  # extract far side point from line (and radius)
            lat_fsp = fsp[0]
            lon_fsp = fsp[1]
            alpha = fsp[2] + 180

            # print(lat_fsp,lon_fsp,alpha)
            for midlati, midlong in zip(midlats,
                                        midlons):  # actual finding corres point (first skip) with helper method get_first_skip
                first_skip = get_first_skip(lat_fsp, lon_fsp, midlati, midlong, alpha / 2)
                first_lat = first_skip
                new_lats.append(first_lat)
                # first_lon = first_skip[1]  # don't include lons
                # new_lons.append(first_lon)

            ff = open(resulting_corresponding_first_skip,
                      'a')  # writing corres points into resulting_corresponding_first_skip
            # for lat, lon in zip(new_lats, new_lons):
            #     ff.write(str(lat) + " " + str(lon) + "\n")

            for lat in zip(new_lats):
                # print(math.degrees(lat[0]))
                ff.write(str(math.degrees(
                    lat[0])) + "\n")  # "lat" and "lon" are tuples, so access the first (and only) element of each
                # convert the lat and lon into degrees, then to string, then write onto file

            ff.write(line)
            ff.write("\n")
            ff.close()

            lats = []
            lons = []  # to reset
            midlats = []
            midlons = []  # to reset

            new_lats = []
            new_lons = []  # to reset

    f.close()
    print("--- end of find_second_skip ---")


def find_third_skip(coors_three, resulting_corresponding_second_skip):
    open(resulting_corresponding_second_skip, 'w').close()  # to clear all present text

    f = open(coors_three, 'r')

    lats = []
    lons = []
    midlats = []
    midlons = []  # originals (reset for each group of points)

    new_lats = []
    new_lons = []  # points we're looking for (reset for each group of points)

    for line in f:
        line = line.strip()

        if "l" not in line:  # sectioning for diff groups of points
            lat, lon, midlat, midlon = line.split(' ')
            lats.append(float(lat))
            lons.append(float(lon))
            midlats.append(float(midlat))
            midlons.append(float(midlon))

        elif "l" in line:  # sectioning for diff groups of points
            fsp = get_coors_from_str(line)  # extract far side point and alpha from line
            lat_fsp = fsp[0]
            lon_fsp = fsp[1]
            alpha = fsp[2] + 180  # Need to +180 since read in 'radius' only

            # print(lat_fsp, lon_fsp, alpha)
            for lati, long, midlati, midlong in zip(lats,
                                                    lons, midlats,
                                                    midlons):  # actual finding corres point (second skip) with helper method get_second_skip

                if midlong < -90:
                    midlong += 360

                # second_skip = get_second_skip(midlati, midlong, lati, long,
                #                               alpha / 2)  # finding antipodes of midpoints, new alpha is half its original
                # second_lat = second_skip
                second_lat = midlati
                new_lats.append(float(second_lat))  # recording the second skip lat

                # second_lon = second_skip[1]
                # new_lons.append(second_lon)  # recording the second skip lon (no need)

            ff = open(resulting_corresponding_second_skip,
                      'a')  # writing corres points into resulting_corresponding_second_skip
            # for lat, lon in zip(new_lats, new_lons):
            #     ff.write(str(lat) + " " + str(lon) + "\n")

            for lat in zip(new_lats):
                # ff.write(str(math.degrees(lat)) + " " +
                #          str(math.degrees(lon)) + "\n")  # "lat" and "lon" are tuples, so access the first (and only) element of each
                # convert the lat and lon into degrees, then to string, then write onto file
                # print(lat, math.degrees(lat[0]))
                ff.write(str(lat[0]) + "\n")

            ff.write(line + "\n")
            ff.close()

            lats = []
            lons = []
            midlats = []
            midlons = []  # to reset

            new_lats = []
            new_lons = []  # to reset

    f.close()
    print("--- end of find-3rd-skip ---")


# helper method (actual math) for find_second_skip; returns list of two numbers: lat, lon
# (https://www.movable-type.co.uk/scripts/latlong-vectors.html)
# lat_fsp = lat far side point
# lon_fsp = lon far side point
# lat_nsp = lat near side point
# lon_nsp = lon near side point
def get_second_skip(lat_fsp, lon_fsp, lat_nsp, lon_nsp, alpha):
    v_1 = np.array([math.cos(math.radians(lon_fsp)) * math.cos(math.radians(lat_fsp)),
                    math.sin(math.radians(lon_fsp)) * math.cos(math.radians(lat_fsp)), math.sin(math.radians(lat_fsp))])
    v_2 = np.array([math.cos(math.radians(lon_nsp)) * math.cos(math.radians(lat_nsp)),
                    math.sin(math.radians(lon_nsp)) * math.cos(math.radians(lat_nsp)), math.sin(math.radians(lat_nsp))])

    c = 1 / (2 * math.cos(math.radians(1 / 3 * alpha)))

    v_3 = numpy.add(v_2, c * v_1)

    v_4 = v_3 / np.linalg.norm(v_3)  # normalizing v_3

    # return [math.asin(v_4[2]), math.asin(v_4[1] / math.cos(math.asin(v_4[2])))
    return math.asin(v_4[2])  # return only desired lat in radian


# actual math of find_first_skip (just halfway "cut")
def get_first_skip(lat_fsp, lon_fsp, lat_nsp, lon_nsp, alpha):
    v_1 = np.array([math.cos(math.radians(lon_fsp)) * math.cos(math.radians(lat_fsp)),
                    math.sin(math.radians(lon_fsp)) * math.cos(math.radians(lat_fsp)), math.sin(math.radians(lat_fsp))])
    v_2 = np.array([math.cos(math.radians(lon_nsp)) * math.cos(math.radians(lat_nsp)),
                    math.sin(math.radians(lon_nsp)) * math.cos(math.radians(lat_nsp)), math.sin(math.radians(lat_nsp))])

    c = 1 / (2 * math.cos(math.radians(1 / 3 * alpha)))

    v_3 = numpy.add(v_2, c * v_1)

    v_4 = v_3 / np.linalg.norm(v_3)  # normalizing v_3

    return math.asin(v_4[2])  # return only desired lat in radian


# helper method for get_second_skip and get_first_skip
# returns the far side point from a string such as "far lat: 0, far lon: 0, radius: 40.10000000000032"
def get_coors_from_str(line):
    split_str = line.split(" ")

    # print(split_str)
    lat = float(split_str[2])  # remove non-digit chars
    lon = float(split_str[5])
    current_radius = float(split_str[7])

    return [lat, lon, current_radius]


def get_mid_corres_normal(lat_fsp, lon_fsp, lat_nsp, lon_nsp):
    v_1 = np.array([math.cos(math.radians(lon_fsp)) * math.cos(math.radians(lat_fsp)),
                    math.sin(math.radians(lon_fsp)) * math.cos(math.radians(lat_fsp)), math.sin(math.radians(lat_fsp))])
    v_2 = np.array([math.cos(math.radians(lon_nsp)) * math.cos(math.radians(lat_nsp)),
                    math.sin(math.radians(lon_nsp)) * math.cos(math.radians(lat_nsp)), math.sin(math.radians(lat_nsp))])

    v_3 = [v_1[0] + v_2[0],
           v_1[1] + v_2[1],
           v_1[2] + v_2[2]]

    v_4 = v_3 / np.linalg.norm(v_3)  # normalizing v_3

    return math.degrees(math.asin(v_4[2]))  # return only desired latitude (but should I returns both?)


# coors contains groups of lat coordinates for different far side point's antipodes
def count_lat_bands(coors, resulting_lat_band_counts):
    open(resulting_lat_band_counts, 'w').close()  # to clear all present text

    ## The lon range for the far side is [-90, +90], with each grid tick bring 0.8,
    ## there are 225 (= 180/0.8) points to calculate for any given lat
    ## Since the calculation for lon is symmetric, we just need to
    ## multiply this number (225) without any further computation.
    # number_long_grid = 225/11.0 # because we did 11 different far-lons; so we need to multiply the same
    number_long_grid = 1

    f = open(coors, 'r')

    bands = [0] * 37
    lats = []  # store lats for each group of points (resets each time)

    ff = open(resulting_lat_band_counts, 'a')  # to record band counts for each far side point's antipode

    for line in f:
        line = line.strip()

        if "l" not in line:  # sectioning for diff groups of points
            if " " in line:
                lats.append(float(line.split(" ")[0]))
            else:
                lats.append(float(line))

        elif "l" in line:  # sectioning for diff groups of points
            far_lat = get_coors_from_str(line)[0]

            for lat in lats:
                lat = float(lat)

                print(str(int(lat / 5 + 18)))
                bands[int(lat / 5 + 18)] += 1  # determining which band this point is in

            ff.write("far lat: " + str(far_lat) + " " + str(bands))
            bands = [0] * 24  # -don't need reset because we add all bands together- (commented back in)
            lats = []  # to reset

    bands = [int(element * number_long_grid) for element in bands]

    ff.write(str(bands))
    print(bands)
    ff.close()


# find midpoint between exhautive-found point and its fsp to record after each point in data1
# find middle with twotwo method cut in half, if atan2 <-90 then add 360, then find antipode (not implemented here)
def get_midpoint(lat_fsp, lon_fsp, lat_nsp, lon_nsp):
    v_1 = np.array([math.cos(math.radians(lon_fsp)) * math.cos(math.radians(lat_fsp)),
                    math.sin(math.radians(lon_fsp)) * math.cos(math.radians(lat_fsp)), math.sin(math.radians(lat_fsp))])
    v_2 = np.array([math.cos(math.radians(lon_nsp)) * math.cos(math.radians(lat_nsp)),
                    math.sin(math.radians(lon_nsp)) * math.cos(math.radians(lat_nsp)), math.sin(math.radians(lat_nsp))])

    v_3 = [v_1[0] + v_2[0],
           v_1[1] + v_2[1],
           v_1[2] + v_2[2]]

    v_4 = v_3 / np.linalg.norm(v_3)  # normalizing v_3

    mid_lat = math.degrees(math.asin(v_4[2]))
    mid_lon = math.degrees(math.atan2(v_4[1], v_4[0]))

    # because the angle as formed by v_1 and v_2 is > 180, we need to find the antipod of the mid-point
    # print([mid_lat, mid_lon])
    mid_lon += 180
    if mid_lon > 270:
        mid_lon -= 360
    mid_lat *= -1
    # print([mid_lat, mid_lon])
    return [mid_lat, mid_lon]
    # return [math.degrees(math.asin(v_4[2])), math.degrees(math.asin(v_4[1] / math.cos(math.asin(v_4[2]))))]



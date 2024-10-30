
"""
This file is part of the OCULARIS Ocular Proton Therapy Treatment Planning System. 
It is subject to the license terms in the LICENSE file located in the top-level directory of this distribution.

This program is not certified for clinical use and is provided WITHOUT ANY WARRANTY or implied warranty.
For accuracy, users should validate OCULARIS independently before drawing any conclusions.

"""


from collections.abc import Iterable
from typing import Tuple
import numpy as np
from numba import jit
from PIL import Image, ImageDraw
from math import isclose, floor
import pyvista as pv
import copy

from pyproton.volume import Grid
from pyproton.structure.organ_data import OrganData
from pyproton.structure.structure import Structure, StructureSlice
from pyproton.utils.concave_hull import ConcaveHull


class ArrayStructure(Structure):
    def __init__(self, organ: OrganData, array: np.ndarray, grid: Grid = None):
        self._organ = organ
        self._array = np.asarray(sorted(array, key=lambda x: x[2]))
        self._points_outside_grid = np.zeros([0,3])
        self._grid = grid
        self._mask = None
        self._points_of_binary_mask = None
        self._color = None
        self._points_resampled = False
        self._mask_defined_already = False
        self._title = f"Structure {self._organ.name} of type {self._organ.type}"

    @property
    def organ(self) -> OrganData:
        return self._organ

    @property
    def color(self) -> np.ndarray:
        if (np.asarray(self._color) == None).any():
            self._color = np.asarray([1, 1, 1])
        return self._color

    @color.setter
    def color(self, color):
        self._color = color

    @property
    def grid(self) -> Grid:
        return self._grid

    @grid.setter
    def grid(self, grid):
        self._grid = grid

    @property
    def com(self):
        #centre of mass
        return np.asarray([sum(row[i] for row in self.contour_coordinates)/len(self.contour_coordinates) for i in range(len(self.contour_coordinates[0]))])

    @property
    def contour_slices(self) -> Iterable[StructureSlice]:
        raise NotImplementedError()

    @property
    def contour_coordinates(self) -> np.ndarray:
        return self._array

    @contour_coordinates.setter
    def contour_coordinates(self, new_coordinates):
        self._array = new_coordinates
        
    @property
    def structure_points(self) -> np.ndarray:
        if len(self._points_outside_grid) != 0:
            return np.concatenate(
                (self._array, self._points_outside_grid))
            
        return self.contour_coordinates     
    

    
    @property
    def binary_mask(self) -> np.ndarray:

        if self._mask is None:
            print(f"Defining binary mask for {self._organ.name}")
            self._mask = _binary_mask_from_contours(
                self._array,
                self._grid.size,
                self._grid.spacing,
                self._grid.origin,
                )
            self._mask = np.swapaxes(self._mask, 0, 1)
            self._mask_defined_already = True
        return self._mask
    
    @binary_mask.setter
    def binary_mask(self, mask):
        self._mask = mask
        print(f"Binary mask set for {self._organ.name}")
            

    @property
    def binary_mask_regenerate(self) -> np.ndarray:
        print(f"Defining binary mask for {self._organ.name}")
        self._mask = _binary_mask_from_contours(
            self._array,
            self._grid.size,
            self._grid.spacing,
            self._grid.origin,
            )
        self._mask = np.swapaxes(self._mask, 0, 1)
        return self._mask
    
    @property
    def points_of_binary_mask(self):
        
        self.resample_contour_points_in_grid(1)
        
        xes, yes, zes = np.where(self.binary_mask)
        
        my_indices = list(zip(xes,yes,zes))
        
        points = np.zeros(shape=(len(my_indices), len(my_indices[0])))
        
        for idx in range(len(my_indices)):
            x_bin, y_bin, z_bin = my_indices[idx]
        
            x_pos = self.grid.spacing[0] * \
                (x_bin + 0.5) + self.grid.origin[0]
            y_pos = self.grid.spacing[1] * \
                (y_bin + 0.5) + self.grid.origin[1]
            z_pos = self.grid.spacing[2] * \
                (z_bin + 0.5) + self.grid.origin[2]
        
            points[idx,0:] = np.asarray([x_pos, y_pos, z_pos])   
        
 
        self._points_of_binary_mask = points
        
        return self._points_of_binary_mask
    

    @property
    def volume(self) -> float:
        volume_of_voxel = self._grid.spacing[0] * \
            self._grid.spacing[1]*self._grid.spacing[2]
        return np.sum(self.binary_mask)*volume_of_voxel

    @property
    def bounding_box(self):
        # TODO: Specify return type and format.
        raise NotImplementedError()

    @property
    def convex_hull(self):
        # TODO: Specify return type and format.
        raise NotImplementedError()

    def __repr__(self):
        return self._title

    def resample_contour_points_in_grid(self, sampling_density=1, force = 0):
        
        if self._points_resampled == True and force == 0:
            # print("Structure already resampled. Did nothing")
            return
        
        
        cloud = pv.PolyData(self.contour_coordinates)
        # cloud.plot()

        volume = cloud.delaunay_3d(alpha=14.)
        shell = volume.extract_geometry()
        # shell.plot()

        points_outside_grid = np.asarray(list(filter(lambda x: x[2] < self.grid.origin[2] or x[2] > self.grid.origin[2] + self.grid.spacing[2]*self.grid.size[2]
                                                     or x[1] < self.grid.origin[1] or x[1] > self.grid.origin[1] + self.grid.spacing[1]*self.grid.size[1]
                                                     or x[0] < self.grid.origin[0] or x[0] > self.grid.origin[0] + self.grid.spacing[0]*self.grid.size[0], self.contour_coordinates)))

        new_points = np.empty((0, 3), dtype=np.float64)
        new_points = np.concatenate(
            (self.contour_coordinates, new_points))
        for z_step in np.linspace(0, self.grid.size[2], int(self.grid.size[2]*sampling_density)):
            z = self.grid.origin[2] + self.grid.spacing[2]*z_step
            resampled_points = _resample_points_in_plane(
                shell, z, self.grid)
            new_points = np.concatenate((new_points, resampled_points))
        
        if len(new_points) == 0:
            print("Warning. Not resampled")
            return
        
        self.contour_coordinates = new_points
        

        if points_outside_grid.shape[0] > 0:
            self._points_outside_grid = points_outside_grid

        self._points_resampled = True
        print(f"Contour points of {self._organ.name} resampled in grid")


def items_in_list_within_tolerance(list_of_values, value, atol):

    return np.isclose(
        np.asarray(list_of_values),
        value,
        atol=atol,
    ).all()


def _resample_points_in_plane(shell, z, grid, atol = 1e-3):
    
    points_in_plane = []

    normal = (0, 0, 1)
    # z_pos = grid.spacing[2] * (z + 0.5) + grid.origin[2]

    clipped = shell.clip(normal=normal, origin=[0, 0, z])
    # clipped.plot()

    cell_indices = clipped.faces

    cell_node_ids = cell_indices.reshape(-1, 4)[:, 1:4].ravel()

    cell_nodes = clipped.points[cell_node_ids]

    filtered_points = list(filter(lambda p: abs(p[2] - z) < atol, cell_nodes))

    for p in filtered_points:
        points_in_plane.append([p[0], p[1], z])

    if len(filtered_points) == 0:
        return np.empty((0, 3))

    return np.array(points_in_plane)


def _binary_mask_from_contours(contours: np.ndarray,
                               grid_size: Tuple[int, int, int],
                               grid_spacing: np.ndarray,
                               grid_origin: np.ndarray) -> np.ndarray:
    '''
    Calculates the binary 3D mask for an entire structure.

    Parameters
    ----------
    contours: np.ndarray
        Contours to convert to a binary mask
    grid_size: tuple of int
        Shape: (3,)
    grid_spacing: np.ndarray
        Shape: (3,)
    grid_origin: np.ndarray
        Shape: (3,)

    Returns
    -------
    np.ndarray
        Binary mask of shape grid_size.
    '''

    floating_point_fudge_factor = 0.99999999

    plane_masks = np.full(grid_size, False)
    contours = _split_slices(contours, grid_spacing[2], grid_origin[2])
    for contour in contours:
        assert items_in_list_within_tolerance(
            contour[:, 2], contour[0, 2], grid_spacing[2])
        z = contour[0, 2]
        z_index = int(floating_point_fudge_factor *
                      (z - grid_origin[2]) // grid_spacing[2])
        if z_index >= grid_size[2]:
            continue
        mask = _contour_slice_to_mask(contour[:, :2],
                                      grid_origin[:2],
                                      grid_spacing[:2],
                                      grid_size[:2])
        plane_masks[:, :, z_index] = mask

    return plane_masks


def _split_slices(points: np.ndarray, slice_spacing: float, z_origin: float):
    '''
    Take an array and split it into a list of sub-arrays based on the
    z-coordinate.

    Parameters
    ----------
    points: np.ndarray
        Array representing a set of points. The points are assumed to represent
        a structure contour, i. e. subsequent points within the same slice are
        supposed to be connected by lines. Within a slice, all points share the
        same z-coordinate. The points are supposed to be ordered such that
        the points belonging to the same slice form a contiguous subarray.
        Shape: (n_points, 3)
    slice_spacing: float
        Spacing between slices. This is needed when there was a registration
        and the points need to be mapped to slices first.
    z_origin: float
        z-component of the coordinate system origin. The slice with index 0
        has this z-coordinate.

    Returns
    -------
    list of np.ndarray
        Each element in the list represents the points in a single slice.
        Empty slices are omitted.
    '''
    slices = []
    z_indices = np.rint((points[:, 2] - z_origin) /
                        slice_spacing).astype(np.int64)

    i_start = 0
    i_end = 0
    while i_end < points.shape[0]:
        if z_indices[i_end] != z_indices[i_start]:
            # Save the old slice and start a new one.
            slices.append(points[i_start:i_end, :])
            i_start = i_end
        i_end += 1

    assert i_end == points.shape[0]

    if i_end - i_start > 1:
        # Add the last slice.
        slices.append(points[i_start:i_end, :])

    return slices


def _unscramble_contour_points(points):

    length_of_ruler = 80

    hull = ConcaveHull()
    hull.loadpoints(points)

    hull.calculatehull(length_of_ruler)
        


    xes, yes = hull.boundary.boundary.xy

    unscrabled_points = []

    for x, y in zip(xes, yes):
        unscrabled_points.append((x, y))

    return unscrabled_points


def _contour_slice_to_mask(points: np.ndarray,
                           origin: np.ndarray,
                           spacing: np.ndarray,
                           size: Tuple[int, int]):
    '''
    Returns the indices of grid cells that are intersected by the given contour
    points.

    This function assumes that the input points lie in the same CT slice,
    i. e. they share the same z coordinate. The z-coordinate must be the
    coordinates of the slice exactly.

    Parameters
    ----------
    points: np.ndarray
        The points defining the contour. Shape: (n_points, 2)
    origin: np.ndarray
        Origin of the grid. Shape: (2,)
    spacing: np.ndarray
        Spacing of the grid. Shape: (2,)
    size: tuple of int
        Number of cells in the grid in each direction. Shape: (2,)

    Returns
    -------
    np.ndarray
        The binary mask of the grid. Cells whose value is True are intersected
        by the contour, the others are not.
        A cell is considered to be intersected by the contour if at least
        half of the 16 subsquares are filled by Pillows polygon fill algorithm.
        Shape: size

    Notice
    ------
    The input points are assumed to be in the correct order, i. e. we assume
    that the contour is generated by connecting subsequent rows, and closing
    the contour in the end (no need to add the first point again at the end).
    '''

    # Increase the resolution of the rasterisation by making the grid finer
    # and then average over the superfluous voxels.
    superresolution_factor = 4
    assert np.isclose((superresolution_factor / 2) % 1, 0)
    size = (size[0] * superresolution_factor, size[1] * superresolution_factor)
    spacing = spacing / superresolution_factor

    points = list(map(tuple, (points - origin) / spacing))

    # Removes duplicates
    points = list(set(points))

    if len(points) == 1:
        binary_array = np.zeros([size[0], size[1]], dtype=bool)
        point = points[0]
        idx = floor(abs(origin[0] - point[0])/spacing[0])
        idy = floor(abs(origin[1] - point[1])/spacing[1])
        if idx < binary_array.shape[0] and idy < binary_array.shape[1]:
            binary_array[idx, idy] = True
        return np.round(_average_superresolution(binary_array, superresolution_factor)).astype(np.bool)

    if len(points) >= 3:
        points = _unscramble_contour_points(points)

    # Implementation stolen from (2021-07-20):
    # https://github.com/KeremTurgutlu/dicom-contour/blob/master/RT2MASK.ipynb
    # 'L': 8-bit pixels, black-and-white, one pixel per byte
    # We initialise the image with color zero, i. e. black.
    img = Image.new('L', (size[1], size[0]), color=0)
    # Draw the polygon, the outline will be drawn in black and the fill will
    # have color 1.

    ImageDraw.Draw(img).polygon(points, outline=1, fill=1)
    mask = np.array(img).astype(bool)

    return np.round(_average_superresolution(mask, superresolution_factor)).astype(np.bool)


@jit(nopython=True)
def _average_superresolution(mask: np.ndarray, superres_factor: float):
    # Assumption: The superresolution factor is even.
    # assert np.isclose((superres_factor / 2) % 1, 0)
    new_mask = np.zeros(
        (mask.shape[0] // superres_factor, mask.shape[1] // superres_factor)
    )
    for i in range(new_mask.shape[0]):
        for j in range(new_mask.shape[1]):
            old_i_range = np.arange(i*superres_factor, (i+1)*superres_factor)
            old_j_range = np.arange(j*superres_factor, (j+1)*superres_factor)

            for old_i in old_i_range:
                for old_j in old_j_range:
                    new_mask[i, j] += mask[old_i, old_j]
            new_mask[i, j] /= (superres_factor * superres_factor)

    return new_mask

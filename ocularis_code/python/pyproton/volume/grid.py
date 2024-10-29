from dataclasses import dataclass
import numpy as np


class MeshGrids():
    """
    Use cases for plotting:
    
    image = dcm_volume.pixel_array[0:,0:, z_bin]
    plt.pcolor(structure1.grid.meshgrids.trans_X, structure1.grid.meshgrids.trans_Y, image, cmap=plt.cm.bone)
    
    image = dcm_volume.pixel_array[0:,y_bin, 0:]
    plt.pcolor(structure1.grid.meshgrids.Z, structure1.grid.meshgrids.X, image, cmap=plt.cm.bone)    
    
    image = dcm_volume.pixel_array[x_bin,0:, 0:]
    plt.pcolor(structure1.grid.meshgrids.Z, structure1.grid.meshgrids.Y, image, cmap=plt.cm.bone)    
    """

    def __init__(self, grid):
        self.grid = grid
        self.trans_X, self.trans_Y = [], []
        self.Z, self.X, self.Y = [], [], []

        self.update(self.grid)

    def update(self, grid):
        self.grid = grid

        self.mesh_apex = np.asarray([self.grid.origin[0] + self.grid.size[0]*self.grid.spacing[0], self.grid.origin[1] +
                                    self.grid.size[1]*self.grid.spacing[1], self.grid.origin[2] + self.grid.size[2]*self.grid.spacing[2]])

        self.Z, self.X = np.meshgrid(np.arange(self.grid.origin[2], self.mesh_apex[2] + self.grid.spacing[2], self.grid.spacing[2]), np.arange(
            self.grid.origin[0], self.mesh_apex[0] + self.grid.spacing[0], self.grid.spacing[0]))

        _, self.Y = np.meshgrid(np.arange(self.grid.origin[2], self.mesh_apex[2] + self.grid.spacing[2], self.grid.spacing[2]), np.arange(
            self.grid.origin[1], self.mesh_apex[1] + self.grid.spacing[1], self.grid.spacing[1]))

        self.trans_X, self.trans_Y = np.meshgrid(np.arange(self.grid.origin[0], self.mesh_apex[0] + self.grid.spacing[0], self.grid.spacing[0]), np.arange(
            self.grid.origin[1], self.mesh_apex[1] + self.grid.spacing[1], self.grid.spacing[1]))


@dataclass
class Grid:
    origin: np.ndarray
    spacing: np.ndarray
    size: np.ndarray
    unit: str = "cm"

    def __post_init__(self):
        self._meshgrids = MeshGrids(self)

    @property
    def meshgrids(self):
        self._meshgrids.update(self)
        return self._meshgrids

    def is_finer_than(self, other) -> bool:
        finer = [self.spacing[i] <= other.spacing[i] for i in range(3)]
        return all(finer)

    def __eq__(self, other) -> bool:
        atol = 1e-6

        if not isinstance(other, Grid):
            raise ValueError

        for i in range(3):
            if self.size[i] != other.size[i]:
                return False
            if abs(self.origin[i] - other.origin[i]) > atol:
                return False
            if abs(self.spacing[i] - other.spacing[i]) > atol:
                return False

        if self.unit != other.unit:
            return False

        return True

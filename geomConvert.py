
import damask

grid = damask.GeomGrid.load_ASCII('grainIDs.geom')

grid.save('grainsID')
#
# Author: John Grime 2021 (Emerging Technologies, U. Oklahoma Libraries)
# Alpha test : not for redistribution without explicit permission!
#
# Script for reading a SHAPEFile definition (i.e., the outline of a shape)
# and extruding it into a 3D slab. The shapes are simplified if needed to
# a specific tolerance by removing points that are "too close together".
# This is also useful for the Triangle library used to triangulate the
# surfaces; this library is a little bit unstable, and can crash when e.g.
# duplicate points exist close to one another. This situation can exist with
# shape data from various places, including local transport departments etc,
# and so some form of simplification (even with very small tolerances) is
# usually enough to remove those issues.
#

import math, sys

#
# We use some non-default modules:
#
#   pip3 install pyshp triangle matplotlib
#
# Matplotlib only used if needed, so imported later on.
#

import shapefile
import triangle as tr

#
# Simple utility class to load th shape data and simplify it.
#

class ShapeUtil:
	def __init__(self, fpath: str = None, min_dr: float = 0):
		self.shapefile, self.shapes = None, None
		if fpath != None: self.load_file(fpath,min_dr)

	def load_file(self, file_path: str, min_dr: float = 0):
		'''
		Loads a shape file, and (potentially) simplifies the outline.
		'''

		self.shapefile = shapefile.Reader(file_path)
		self.shapes = self.shapefile.shapes()

		print()
		print(f'{len(self.shapes)} shapes total. Fields:')
		for f in self.shapefile.fields: print('  ', f)
		print()

		# Simplify?
		if min_dr > 0:
			print(f'Simplifying shapes with min_dr = {min_dr}:')
			for i,shape in enumerate(self.shapes):
				parts, points = self.simplify_shape(i, min_dr)
				a, b = len(shape.parts), len(shape.points)
				c, d = len(parts), len(points)
				print(f'  {i:5d} : {a:5d} {b:6d} => {c:5d} {d:6d}')
				shape.parts, shape.points = parts, points
			print()

		# Skip initial deletion flag field; it doesn't show up in the results
		# of calling self.shapefile.record(i), so breaks get_shape_indices()!
		self.field_keys = [f[0] for f in self.shapefile.fields[1:]]

	def simplify_shape(self, shape_i: int, min_dr: float):
		'''
		Simplifies a shape by ensuring that consecutive points in the shape
		definition are separated by at least min_dr.
		'''

		shape = self.shapes[shape_i]
		parts, points = shape.parts, shape.points
		n_parts, n_points = len(parts), len(points)

		min_dr2 = min_dr*min_dr
		parts_, points_ = [], []

		# Simplify each part of the shape independently.
		for part_i in range(n_parts):

			# Start/stop indices of points representing this part of the shape
			i = parts[part_i]
			j = parts[part_i+1] if part_i<n_parts-1 else n_points

			p0 = points[i]
			parts_.append( len(points_) )
			points_.append(p0) # AFTER update of parts_

			# Note: final point in each part is just the original point again,
			# so skip final point and explicitly close polygon after the loop.
			for p in points[i+1:j-1]:
				# If current point is sufficient distance from last point ...
				dx, dy = p[0]-p0[0], p[1]-p0[1]
				if (dx*dx + dy*dy) < min_dr2: continue
				points_.append(p)
				p0 = p

			# Explicitly close part by reconnecting to the first simplified
			# point in the current part.
			i = parts_[-1]
			points_.append( (points_[i][0],points_[i][1]) )

		return parts_, points_

	def get_shape_indices(self, key: str, vals: [str]):
		'''
		Returns the indices of all shapes with one of specified values under
		record entry "key"
		'''

		shape_indices = []

		if (key == None) or (vals == None):
			print(f'ShapeUtil.get_fields(): no key or value array specified!')
			return shape_indices

		try:
			key_idx = self.field_keys.index(key)
		except ValueError:
			print(f'Key "{key}" does not exist in shape data; possible keys are: {self.field_keys}')
			return shape_indices

		for i,shape in enumerate(self.shapes):
			rec = self.shapefile.record(i) # see note re. field_keys in load_file()
			if rec[key_idx] in vals: shape_indices.append(i)

		return shape_indices


def extrude(vtx, tri, exterior_vtx_idx, seal=True):
	'''
	Return vertex and triangle data required to extrude a planar mesh with
	the specified peripheral points into 3d wedge.
	Note: no z coords are added to the vertices! The data returned is therefore
	still planar until the user adds z coords to the vertices!
	'''

	N, N_idx = len(vtx), len(exterior_vtx_idx)

	#
	# Duplicate all vertices, flip vertex index order for new triangles to
	# "mirror" existing planar surface - ready to add z offset later.
	#

	vtx_ = [ [x,y] for x,y in vtx ]
	tri_ = [ [i+N,j+N,k+N] for i,j,k in tri ]

	#
	# Connect "edges" of old and new planar surface
	#
	# d--c => "old" points
	#    |
	# a--b => "new" points
	#
	
	for i in range(N_idx-1):
		a,b = N+i, N+(i+1)
		d,c = exterior_vtx_idx[i], exterior_vtx_idx[i+1]
		tri_.extend( [ [a,b,c], [c,d,a] ] )

	# If sealing the extrusion, connect "final" points to "first" points.
	if seal == True:
		a,b = N+(N_idx-1), N
		d,c = exterior_vtx_idx[-1], exterior_vtx_idx[0]
		tri_.extend( [ [a,b,c], [c,d,a] ] )

	return vtx_, tri_


#
# Simplify shape data via min point separation.
# Also used to remove adjacent duplicate points, which can crash
# the "triangle" library!
#

county_key = 'COUNTY_NAM'
min_dr = 1e-6
plot_shapes = False

#
# Get user parameters
#

if len(sys.argv) < 2:
	print('')
	print('Usage:')
	print(f'  python3 {sys.argv[0]} SHAPEFile prefix [shape name] [delta z]')
	print('')
	print('Where:')
	print('  SHAPEFile  : path prefix to SHAPEFile data (don\'t include the file suffix!)')
	print('  shape name : name of shape in file (if unspecified, shapes in file listen then exit)')
	print('  delta z    : extrusion height on z axis (default: 1.0)')
	print('')
	print('Example:')
	print(f'  python3 {sys.argv[0]} my_data/my_shapefile blah 0.25')
	print('')
	print('  This looks for the "blah" shape in the my_data/my_shapefile.(dbf,shp,shx) file')
	print('  combination, and extrudes it into a 3d model of "height" 0.25 units. The output')
	print('  is written to an "output.obj" file.')
	print('')
	sys.exit(-1)

fprefix = sys.argv[1]
shape_name = sys.argv[2] if (len(sys.argv)>2) else None
dz = float(sys.argv[3]) if (len(sys.argv)>3) else 1.0

#
# Load shape data from file
#

su = ShapeUtil(fprefix,min_dr=min_dr)

#
# If no target county name specified, print the shape data and exit.
#

if shape_name == None:
	# check we have the appropriate county name field
	try:
		_ = su.field_keys.index(county_key)
	except:
		print(f'Unable to find data under field key "{county_key}" for input file!')
		sys.exit(-1)

	for i,s in enumerate(su.shapes):
		rec = su.shapefile.record(i)
		print( ' '.join([f'{k}={rec[k]}' for k in su.field_keys]) )
	
	sys.exit(0)

#
# Get specific county data
#

shape_indices = su.get_shape_indices(key=county_key,vals=[shape_name])
if len(shape_indices) < 1:
	print(f'Unable to find shape "{shape_name}"')
	sys.exit(-1)

shape_idx = shape_indices[0]
shape = su.shapes[shape_idx]
parts, points = shape.parts, shape.points

print(f'{shape_name} ({shape_idx}): {len(parts)} parts, {len(points)} points')

#
# Break up longer line sections so they do not exceed 2 x min_dr? This is now
# redundant, instead we specify triangle quality params for "triangle" library
#

if False:
	N = len(points)
	vtx = [ [x,y] for x,y in points[:1] ]
	max_dr = min_dr*2

	for i in range(1,N):
		v0, v1 = points[i-1], points[i]
		[dx, dy], dz = [v1[j]-v0[j] for j in (0,1)], 0
		dr = math.sqrt(dx*dx + dy*dy + dz*dz)
		if dr > max_dr:
			n = int( math.floor(dr/max_dr) )
			for j in range(n):
				vtx.append( [v0[0]+(dx/n)*(j+1),v0[1]+(dy/n)*(j+1)] )
		vtx.append(v1)

	points = vtx

#
# While the triangle library can add new points etc depending on the arguments,
# it does seems to preserve the initial vertex data order. This means we should
# be able to assume the *initial* N points will all remain on the boundary.
#

# Shapes can "close" via copyig 1st point as final point; remove if so
if points[0] == points[-1]: points = points[:-1]

# All N initial points are assumed to define the shape boundary.
N = len(points)
vtx, tri = points, []

if True:
	# Generate "segments" (i.e., edges); include final vtx -> first vtx
	segments = [[i,i+1] for i in range(N-1)]
	segments.append( [N-1,0] )

	A = dict(vertices=points, segments=segments)
#	B = tr.triangulate(A, 'p' )  # p = no new vertices, just use existing
	B = tr.triangulate(A, 'pq' ) # pq = quality control on angles; "better" mesh

	# Plot results?	
	if plot_shapes == True:
		import matplotlib.pyplot as plt
		tr.compare(plt, A, B)
		plt.show()

	vtx, tri = B['vertices'], B['triangles']

Nv1, Nv2, Nt = len(points), len(vtx), len(tri)
print(f'{Nv1} input vtx; {Nv2} output vtx; {Nt} output tri')

vtx_, tri_ = extrude(vtx, tri, [i for i in range(N)])

# Note: vertices are all 2D at this stage, so add z coords
all_vtx = [[x,y,-dz/2] for x,y in vtx] + [[x,y,+dz/2] for x,y in vtx_]
all_tri = [[k,j,i] for i,j,k in tri] + [[i,j,k] for i,j,k in tri_]

# Write simple Wavefront .obj file
with open('output.obj','w') as f:
	for v in all_vtx:  print(f'v {v[0]} {v[1]} {v[2]}', file=f)	
	for t in all_tri:  print(f'f {t[0]+1} {t[1]+1} {t[2]+1}', file=f)

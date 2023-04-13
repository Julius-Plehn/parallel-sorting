project('parallel-sorting', 'c')

cc = meson.get_compiler('c')
mpic = dependency('mpi', language : 'c', method :'config-tool')
omp_dep = dependency('openmp')

deps = [mpic,omp_dep]

sorting_incs = include_directories([
	'include',
])

sorting_srcs = files([
    'lib/sorting.c'
])

parallel_sorting = executable('parallel-sorting', sorting_srcs,
	dependencies: [deps],
	include_directories: [sorting_incs],
)


project('seticore', ['cpp', 'c', 'cuda'],
        default_options: [
            'buildtype=release',
            'cpp_std=c++14',
            'werror=true',
        ])

cuda_args = [ '--std=c++14' ]

cmake = import('cmake')

fmt_subproj = cmake.subproject('fmt')
fmt_dep = fmt_subproj.dependency('fmt')

boost_dep = dependency('boost', modules: [
    'filesystem',
    'program_options',
    'system',
])

hdf5_dep = dependency('hdf5', language: 'c')

cuda_dep = dependency('cuda', version: '>=11', modules: ['cufft'])

capnp_opt = cmake.subproject_options()
capnp_opt.set_override_option('warning_level', '0')
capnp_subproj = cmake.subproject('capnproto', options: capnp_opt)
kj_dep = capnp_subproj.dependency('kj')
capnp_dep = capnp_subproj.dependency('capnp')

deps = [fmt_dep, boost_dep, hdf5_dep, cuda_dep, kj_dep, capnp_dep]

srcs = [
    'beamformer.cu',
    'beamforming_config.cpp',
    'cuda_util.cu',
    'dedoppler.cu',
    'dedoppler_hit.cpp',
    'dat_file_writer.cpp',
    'filterbank_buffer.cu',
    'filterbank_file.cpp',
    'fil_file.cpp',
    'h5_file.cpp',
    'hit.capnp.c++',
    'hit_file_writer.cpp',
    'hit_recorder.cpp',
    'multibeam_buffer.cu',
    'raw_file_group.cpp',
    'recipe_file.cpp',
    'run_dedoppler.cpp',
    'util.cpp',
    'vector_hit_recorder.cpp',
]

tests = [
    'dedoppler_test.cpp',
    'multibeam_buffer_test.cpp',
]

executable('seticore', srcs + ['main.cpp'],
           dependencies: deps, cuda_args: cuda_args)

executable('beamforming_integration_test',
           srcs + ['beamforming_integration_test.cpp'],
           dependencies: deps, cuda_args: cuda_args)

executable('tests', srcs + tests + ['tests.cpp'],
           dependencies: deps, cuda_args: cuda_args)

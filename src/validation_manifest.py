"""Required evidence IDs for the implemented baseline and reserved full profile."""

# This ordered registry freezes the historical check-to-ID correspondence.
LEGACY_SYMBOLIC_LABELS = (
    'three-dimensional quadrupole trace',
    'three-dimensional octupole traces',
    'spin-density inverse map on all 16 Hermitian basis matrices',
    'octupole longitudinal eigenvalues',
    'octupole helicity-one matrix',
    'quadrupole helicity-two matrix',
    'STF momentum tensors through rank five',
    'STF of rank-two product vanishes',
    'dual STF of rank-two product vanishes',
    'dual STF of rank-three double contraction vanishes',
    'quark complex low branch m=1',
    'quark complex high branch m=1',
    'gluon complex low branch m=1',
    'gluon complex high branch m=1',
    'sine invariant m=1',
    'quark complex low branch m=2',
    'quark complex high branch m=2',
    'gluon complex low branch m=2',
    'gluon complex high branch m=2',
    'sine invariant m=2',
    'quark complex low branch m=3',
    'quark complex high branch m=3',
    'gluon complex low branch m=3',
    'gluon complex high branch m=3',
    'sine invariant m=3',
    'quark catalogue has 32 entries',
    'quark 64-by-32 tensor map has rank 32',
    'quark multipole counts 2,6,10,14',
    'gluon catalogue has 32 entries',
    'gluon 64-by-32 tensor map has rank 32',
    'gluon multipole counts 2,6,10,14',
    'quark projection Gram matrix',
    'fragmentation projection Gram matrix',
    'electromagnetic trace channel 0,0',
    'electromagnetic trace channel 0,2',
    'electromagnetic trace channel 0,3',
    'electromagnetic trace channel 1,0',
    'electromagnetic trace channel 1,2',
    'electromagnetic trace channel 1,3',
    'electromagnetic trace channel 2,0',
    'electromagnetic trace channel 2,2',
    'electromagnetic trace channel 2,3',
    'electromagnetic trace channel 3,0',
    'electromagnetic trace channel 3,2',
    'electromagnetic trace channel 3,3',
    'SIDIS scalar-helicity-Collins master contraction',
    'ideal transition vector moment',
    'ideal transition quadrupole moment',
    'ideal transition octupole moment',
 )
LEGACY_IDS = {label: f"legacy.symbolic.{i:03d}" for i, label in enumerate(LEGACY_SYMBOLIC_LABELS, 1)}
BORN_CHECKS = (
    'on_shell_and_conservation', 'photon_Ward', 'gluon_Ward',
    'Hermitian_hard_matrix', 'positive_semidefinite_hard_matrix',
)
SOFTWARE_TESTS = ('test_symbolic_report', 'test_born_report', 'test_grid')
BASELINE_REQUIRED = tuple(LEGACY_IDS.values()) + tuple(f'born.{name}' for name in BORN_CHECKS) + ('grid.complete',) + tuple(f'software.{name}' for name in SOFTWARE_TESTS)

# Suite IDs reserve the scope of the complete manuscript validation. They are
# not synthetic PASS checks. Full remains incomplete until these suites exist.
FULL_PENDING_SUITES = (
    'suite.spin_independent', 'suite.stf_independent', 'suite.gluon_dictionary',
    'suite.correlators_projectors', 'suite.sidis_complete', 'suite.dy_complete',
    'suite.universality', 'suite.gluon_response', 'suite.born_independent_scan',
    'suite.positivity_collinear_fourier', 'suite.negative_controls',
)
PROFILE_REQUIRED = {
    'baseline': BASELINE_REQUIRED,
    'full': BASELINE_REQUIRED + FULL_PENDING_SUITES,
}

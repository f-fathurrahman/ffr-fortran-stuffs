include ../make.inc

LIBNAME = lib78_effpot.a

SRC = \
compute_anharmonics.F90   m_effpot_mpi.F90            m_polynomial_coeff.F90 \
effpot_xml.c              m_anharmonics_terms.F90         m_fit_polynomial_coeff.F90  m_polynomial_conf.F90 \
init10.F90                m_effective_potential.F90       m_harmonics_terms.F90       m_polynomial_term.F90 \
interfaces_78_effpot.F90  m_effective_potential_file.F90  m_multibinit_dataset.F90

include ../make.inc.post


compute_anharmonics.o : compute_anharmonics.F90 m_multibinit_dataset.o m_effective_potential_file.o m_effective_potential.o m_anharmonics_terms.o 
init10.o : init10.F90 
interfaces_78_effpot.o : interfaces_78_effpot.F90 m_effective_potential.o m_multibinit_dataset.o 
m_anharmonics_terms.o : m_anharmonics_terms.F90 m_polynomial_coeff.o 
m_effective_potential.o : m_effective_potential.F90 m_multibinit_dataset.o m_polynomial_coeff.o m_effpot_mpi.o m_harmonics_terms.o m_anharmonics_terms.o m_polynomial_conf.o 
m_effective_potential_file.o : m_effective_potential_file.F90 m_polynomial_term.o m_polynomial_coeff.o m_multibinit_dataset.o m_effective_potential.o m_anharmonics_terms.o m_harmonics_terms.o 
m_effpot_mpi.o : m_effpot_mpi.F90 
m_fit_polynomial_coeff.o : m_fit_polynomial_coeff.F90 m_effective_potential.o m_polynomial_coeff.o 
m_harmonics_terms.o : m_harmonics_terms.F90 
m_multibinit_dataset.o : m_multibinit_dataset.F90 
m_polynomial_coeff.o : m_polynomial_coeff.F90 m_polynomial_term.o 
m_polynomial_conf.o : m_polynomial_conf.F90 
m_polynomial_term.o : m_polynomial_term.F90 

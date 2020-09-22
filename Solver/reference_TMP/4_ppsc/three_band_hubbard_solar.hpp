
#ifndef _PPSC_HAMILTONIANS_THREE_BAND_HUBBARD_SOLAR_HPP
#define _PPSC_HAMILTONIANS_THREE_BAND_HUBBARD_SOLAR_HPP

// -----------------------------------------------------------------------
//
// Three band Hubbard Hamiltonian builder
//
// Author: P. Werner, philipp.werner@gmail.com (2016)
//
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
#include "ppsc/ppsc.hpp"
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
namespace ppsc {
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
namespace hamiltonians {
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
enum interaction_type { ising, kanamori };

// -----------------------------------------------------------------------
template<class HILB, interaction_type INTERACTION>
class three_band_hubbard_solar {

public:

  const interaction_type interaction = INTERACTION;

  // ---------------------------------------------------------------------
  three_band_hubbard_solar(int nt, HILB & hilbert_space) :
    site(0),
    nt(nt),
    mu(0.0),
    // -- parameter vectors
    U(nt+2, 0.0), J(nt+2, 0.0), Bz(nt+2, 0.0),
    E1(nt+2, 0.0), E2(nt+2, 0.0), E3(nt+2, 0.0),
    // -- expectation value vectors
    Q_exp(nt+2), Eint_exp(nt+2),
    n1_exp(nt+2), n2_exp(nt+2), n3_exp(nt+2), n_exp(nt+2),
    m1_exp(nt+2), m2_exp(nt+2), m3_exp(nt+2), m_exp(nt+2),
    d1_exp(nt+2), d2_exp(nt+2), d3_exp(nt+2),
    densdens_intraband_exp(nt+2),
    densdens_interband_exp(nt+2),
    densdens_interband_equalspin_exp(nt+2),
    spin_flip_exp(nt+2), pair_hop_exp(nt+2),
    // -- rotation matrices
    rotN2(15,15),rotdagN2(15,15),
    rotN3(20,20),rotdagN3(20,20),
    rotN4(15,15),rotdagN4(15,15),
    //
    GSN1_exp(nt+2),
    GSN2_exp(nt+2),EXCN2_1_exp(nt+2),EXCN2_2_exp(nt+2),
    GSN3_exp(nt+2),EXCN3_1_exp(nt+2),EXCN3_2_exp(nt+2),
    GSN4_exp(nt+2),EXCN4_1_exp(nt+2),EXCN4_2_exp(nt+2),
    GSN5_exp(nt+2),
    //
    hil_(hilbert_space) {

    // -- Construct basic operators
    c1ua = hil_.c_op_[hil_.flavor(0, 0, 0)];
    c1uc = hil_.c_op_[hil_.flavor(0, 0, 1)];
    c1da = hil_.c_op_[hil_.flavor(0, 1, 0)];
    c1dc = hil_.c_op_[hil_.flavor(0, 1, 1)];

    c2ua = hil_.c_op_[hil_.flavor(1, 0, 0)];
    c2uc = hil_.c_op_[hil_.flavor(1, 0, 1)];
    c2da = hil_.c_op_[hil_.flavor(1, 1, 0)];
    c2dc = hil_.c_op_[hil_.flavor(1, 1, 1)];

    c3ua = hil_.c_op_[hil_.flavor(2, 0, 0)];
    c3uc = hil_.c_op_[hil_.flavor(2, 0, 1)];
    c3da = hil_.c_op_[hil_.flavor(2, 1, 0)];
    c3dc = hil_.c_op_[hil_.flavor(2, 1, 1)];
    //
    n1u = c1uc * c1ua; n1d = c1dc * c1da;
    n2u = c2uc * c2ua; n2d = c2dc * c2da;
    n3u = c3uc * c3ua; n3d = c3dc * c3da;
    //
    d1 = n1u * n1d;
    d2 = n2u * n2d;
    d3 = n3u * n3d;
    //
    n1 = n1u + n1d;
    n2 = n2u + n2d;
    n3 = n3u + n3d;
    n = n1 + n2 + n3;
    //
    m1 = n1u - n1d;
    m2 = n2u - n2d;
    m3 = n3u - n3d;
    m = m1 + m2 + m3;
    //
    densdens_intraband = d1  + d2 + d3;
    densdens_interband = n1 * n2 + n1 * n3 + n2 * n3;
    densdens_interband_equalspin =
        n1u * n2u + n1d * n2d
      + n1u * n3u + n1d * n3d
      + n2u * n3u + n2d * n3d;
    //
    if(interaction == interaction_type::kanamori)
    {
      pair_hop =
           c2uc * c2dc * c1ua * c1da
         + c1dc * c1uc * c2da * c2ua
         + c3uc * c3dc * c2ua * c2da
         + c2dc * c2uc * c3da * c3ua
         + c1uc * c1dc * c3ua * c3da
         + c3dc * c3uc * c1da * c1ua;
      spin_flip =
           c1dc * c2uc * c2da * c1ua
         + c1uc * c2dc * c2ua * c1da
         + c2dc * c3uc * c3da * c2ua
         + c2uc * c3dc * c3ua * c2da
         + c3dc * c1uc * c1da * c3ua
         + c3uc * c1dc * c1ua * c3da;
    }
    // -- Construct rotations
    rotN2.setZero(15,15);
    rotN2.row(0) <<0.00000000,	 0.00000000,	 0.70710678,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.70710678,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000;
    rotN2.row(1) <<0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.70710678,	 0.00000000,	 0.00000000,	 0.00000000,	-0.70710678,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000;
    rotN2.row(2) <<0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	-1.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000;
    rotN2.row(3) <<0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.62034084,	-0.77856531,	 0.00000000,	 0.00000000,	-0.09493843;
    rotN2.row(4) <<0.57735027,	-0.40824829,	 0.00000000,	-0.09681950,	 0.00000000,	 0.70044699,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000;
    rotN2.row(5) <<0.00000000,	 0.00000000,	 0.00000000,	-0.70044699,	 0.00000000,	-0.09681950,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.52913878,	 0.38978795,	 0.00000000,	 0.00000000,	 0.26091668;
    rotN2.row(6) <<0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 1.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000;
    rotN2.row(7) <<0.57735027,	-0.40824829,	 0.00000000,	 0.09681950,	 0.00000000,	-0.70044699,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000;
    rotN2.row(8) <<0.57735027,	 0.81649658,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000;
    rotN2.row(9) <<0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.23495011,	 0.29994462,	 0.00000000,	 0.00000000,	-0.92457107;
    rotN2.row(10)<<0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	-0.70710678,	 0.00000000,	 0.00000000,	 0.00000000,	-0.70710678,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000;
    rotN2.row(11)<<0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 1.00000000,	 0.00000000,	 0.00000000;
    rotN2.row(12)<<0.00000000,	 0.00000000,	 0.00000000,	 0.70044699,	 0.00000000,	 0.09681950,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.52913878,	 0.38978795,	 0.00000000,	 0.00000000,	 0.26091668;
    rotN2.row(13)<<0.00000000,	 0.00000000,	-0.70710678,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.70710678,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000;
    rotN2.row(14)<<0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 1.00000000,	 0.00000000;
    rotdagN2 = rotN2.transpose();
    //
    rotN3.setZero(20,20);
    rotN3.row(0) <<0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 1.00000000,	 0.00000000,	 0.00000000,	 0.00000000;
    rotN3.row(1) <<0.70710678,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	-0.70710678,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000;
    rotN3.row(2) <<0.00000000,	 0.00000000,	 0.70710678,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	-0.70710678,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000;
    rotN3.row(3) <<0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.64409803,	 0.00000000,	 0.00000000,	-0.40824829,	 0.00000000,	 0.00000000,	 0.00000000,	-0.29178370,	 0.00000000,	 0.00000000,	 0.57735027,	 0.00000000,	 0.00000000;
    rotN3.row(4) <<0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.70710678,	 0.00000000,	 0.00000000,	-0.29178370,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	-0.64409803,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000;
    rotN3.row(5) <<0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.81649658,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.57735027,	 0.00000000,	 0.00000000;
    rotN3.row(6) <<0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.70710678,	 0.00000000,	 0.00000000,	 0.29178370,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.64409803,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000;
    rotN3.row(7) <<0.00000000,	 0.00000000,	 0.70710678,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.70710678,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000;
    rotN3.row(8) <<0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	-0.64409803,	 0.00000000,	 0.00000000,	-0.40824829,	 0.00000000,	 0.00000000,	 0.00000000,	 0.29178370,	 0.00000000,	 0.00000000,	 0.57735027,	 0.00000000,	 0.00000000;
    rotN3.row(9) <<0.70710678,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.70710678,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000;
    rotN3.row(10)<<0.00000000,	 0.70710678,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	-0.70710678,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000;
    rotN3.row(11)<<0.00000000,	 0.00000000,	 0.00000000,	-0.70710678,	 0.00000000,	 0.00000000,	-0.29178370,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.64409803,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000;
    rotN3.row(12)<<0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.64409803,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	-0.40824829,	 0.00000000,	 0.00000000,	 0.29178370,	 0.00000000,	 0.00000000,	 0.57735027,	 0.00000000;
    rotN3.row(13)<<0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	-0.70710678,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.70710678,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000;
    rotN3.row(14)<<0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	-0.64409803,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	-0.40824829,	 0.00000000,	 0.00000000,	-0.29178370,	 0.00000000,	 0.00000000,	 0.57735027,	 0.00000000;
    rotN3.row(15)<<0.00000000,	 0.00000000,	 0.00000000,	-0.70710678,	 0.00000000,	 0.00000000,	 0.29178370,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	-0.64409803,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000;
    rotN3.row(16)<<0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.81649658,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.57735027,	 0.00000000;
    rotN3.row(17)<<0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	-0.70710678,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	-0.70710678,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000;
    rotN3.row(18)<<0.00000000,	 0.70710678,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.70710678,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000;
    rotN3.row(19)<<0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 1.00000000;
    rotdagN3 = rotN3.transpose();
    //
    rotN4.setZero(15,15);
    rotN4.row(0) <<  0.00000000,  0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 1.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000;
    rotN4.row(1) <<  0.00000000,  0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 1.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000;
    rotN4.row(2) <<  0.00000000,  0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 1.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000;
    rotN4.row(3) << -0.57735027,  0.00000000,	-0.40824829,	 0.69879657,	-0.02540470,	 0.10506169,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000;
    rotN4.row(4) <<  0.00000000,  0.00000000,	 0.00000000,	 0.10808958,	 0.16424080,	-0.67922132,	 0.00000000,	 0.00000000,	 0.00000000,	 0.16064043,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	-0.68861793;
    rotN4.row(5) <<  0.00000000,  0.00000000,	 0.00000000,	 0.00000000,	 0.68729874,	 0.16619398,	 0.00000000,	 0.00000000,	 0.00000000,	-0.68861793,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	-0.16064043;
    rotN4.row(6) <<  0.00000000,  0.00000000,	 0.00000000,	-0.10808958,	-0.16424080,	 0.67922132,	 0.00000000,	 0.00000000,	 0.00000000,	 0.16064043,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	-0.68861793;
    rotN4.row(7) << -0.57735027,  0.00000000,	-0.40824829,	-0.69879657,	 0.02540470,	-0.10506169,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000;
    rotN4.row(8) <<  0.00000000,  0.70710678,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.70710678,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000;
    rotN4.row(9) <<  0.00000000,  0.00000000,	 0.00000000,	 0.00000000,	-0.68729874,	-0.16619398,	 0.00000000,	 0.00000000,	 0.00000000,	-0.68861793,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	-0.16064043;
    rotN4.row(10)<<  0.00000000, -0.70710678,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.70710678,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000;
    rotN4.row(11)<< -0.57735027,  0.00000000,	 0.81649658,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000;
    rotN4.row(12)<<  0.00000000,  0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 1.00000000,	 0.00000000,	 0.00000000,	 0.00000000;
    rotN4.row(13)<<  0.00000000,  0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 1.00000000,	 0.00000000,	 0.00000000;
    rotN4.row(14)<<  0.00000000,  0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 0.00000000,	 1.00000000,	 0.00000000;
    rotdagN4 = rotN4.transpose();
    //
    Q = operator_type::Identity(hil_);
  }

  // ---------------------------------------------------------------------
  void get_hamiltonian(int tstp, operator_type & Htemp) {

    double u = U[tstp+1];
    double j = J[tstp+1];
    double bz = Bz[tstp+1];

    double e1 = E1[tstp+1];
    double e2 = E2[tstp+1];
    double e3 = E3[tstp+1];

    double half_filling_shift = 0.5*(5*u - 10*j);

    Htemp.InitDiagZero(hil_.ssdim_);

    Htemp = Htemp
      + bz * m
      - mu * n
      + e1 * n1 + e2 * n2 + e3 * n3
      + (u      ) * densdens_intraband
      + (u - 2*j) * densdens_interband
      + (     -j) * densdens_interband_equalspin;

    if(interaction == interaction_type::kanamori)
    {
      Htemp = Htemp - j * (spin_flip + pair_hop);
    }

  }

  // ---------------------------------------------------------------------
  void update_exp_vals(int tstp, operators_type & rho) {
    //
    Q_exp[tstp+1]    = expectation_value(tstp, rho, Q   ).real();
    operator_type Ht; get_hamiltonian(tstp, Ht);
    Eint_exp[tstp+1] = expectation_value(tstp, rho, Ht).real();
    //
    n_exp[tstp+1]    = expectation_value(tstp, rho, n   ).real();
    n1_exp[tstp+1]   = expectation_value(tstp, rho, n1  ).real();
    n2_exp[tstp+1]   = expectation_value(tstp, rho, n2  ).real();
    n3_exp[tstp+1]   = expectation_value(tstp, rho, n3  ).real();
    //
    m_exp[tstp+1]    = expectation_value(tstp, rho, m   ).real();
    m1_exp[tstp+1]   = expectation_value(tstp, rho, m1  ).real();
    m2_exp[tstp+1]   = expectation_value(tstp, rho, m2  ).real();
    m3_exp[tstp+1]   = expectation_value(tstp, rho, m3  ).real();
    //
    d1_exp[tstp+1]   = expectation_value(tstp, rho, d1  ).real();
    d2_exp[tstp+1]   = expectation_value(tstp, rho, d2  ).real();
    d3_exp[tstp+1]   = expectation_value(tstp, rho, d3  ).real();
    //
    densdens_intraband_exp[tstp+1] = expectation_value(tstp, rho, densdens_intraband).real();
    densdens_interband_exp[tstp+1] = expectation_value(tstp, rho, densdens_interband).real();
    densdens_interband_equalspin_exp[tstp+1] = expectation_value(tstp, rho, densdens_interband_equalspin).real();
    spin_flip_exp[tstp+1]   = expectation_value(tstp, rho, spin_flip  ).real();
    pair_hop_exp[tstp+1]    = expectation_value(tstp, rho, pair_hop  ).real();
    //
    // Original basis
    Eigen::MatrixXcd rhoN0(1,1)  ;rhoN0.setZero(1,1)  ;
    Eigen::MatrixXcd rhoN1(6,6)  ;rhoN1.setZero(6,6)  ;
    Eigen::MatrixXcd rhoN2(15,15);rhoN2.setZero(15,15);
    Eigen::MatrixXcd rhoN3(20,20);rhoN3.setZero(20,20);
    Eigen::MatrixXcd rhoN4(15,15);rhoN4.setZero(15,15);
    Eigen::MatrixXcd rhoN5(6,6)  ;rhoN5.setZero(6,6)  ;
    Eigen::MatrixXcd rhoN6(1,1)  ;rhoN6.setZero(1,1)  ;
    //
    // Rotated basis
    Eigen::MatrixXcd rhoN2tilde(15,15);rhoN2tilde.setZero(15,15);
    Eigen::MatrixXcd rhoN3tilde(20,20);rhoN3tilde.setZero(20,20);
    Eigen::MatrixXcd rhoN4tilde(15,15);rhoN4tilde.setZero(15,15);
    //
    // Eigenvalues
    Eigen::VectorXcd eigN2(15);eigN2.setZero(15);
    Eigen::VectorXcd eigN3(20);eigN3.setZero(20);
    Eigen::VectorXcd eigN4(15);eigN4.setZero(15);
    //
    // format
    Eigen::IOFormat myFmt(2, 0, ", ", ";\n", "", "", "[", "]");

    if(tstp>=0)
    {
      //
      // N=0 sector
      rhoN0.block(0,0,1,1)  =rho[tstp].M_[0] ;
      //
      // N=1 sector
      rhoN1.block(0,0,3,3)  =rho[tstp].M_[1] ;
      rhoN1.block(3,3,3,3)  =rho[tstp].M_[4] ;
      GSN1_exp[tstp+1] = rhoN1.trace();
      //
      // N=2 sector
      rhoN2.block(0,0,3,3)  =rho[tstp].M_[2] ;
      rhoN2.block(3,3,9,9)  =rho[tstp].M_[5] ;
      rhoN2.block(12,12,3,3)=rho[tstp].M_[8] ;
      rhoN2tilde=rotdagN2*rhoN2*rotN2;
      GSN2_exp[tstp+1]    = rhoN2tilde.block(6,6,9,9).trace();
      EXCN2_1_exp[tstp+1] = rhoN2tilde.block(1,1,5,5).trace();
      EXCN2_2_exp[tstp+1] = rhoN2tilde.block(0,0,1,1).trace();
      Eigen::ComplexEigenSolver<Eigen::MatrixXcd> eS2(rhoN2,false);
      eigN2 = eS2.eigenvalues();
      //
      // N=3 sector
      rhoN3.block(0,0,1,1)  =rho[tstp].M_[3] ;
      rhoN3.block(1,1,9,9)  =rho[tstp].M_[6] ;
      rhoN3.block(10,10,9,9)=rho[tstp].M_[9] ;
      rhoN3.block(19,19,1,1)=rho[tstp].M_[12];
      rhoN3tilde=rotdagN3*rhoN3*rotN3;
      GSN3_exp[tstp+1]    = rhoN3tilde.block(16,16,4,4).trace();
      EXCN3_1_exp[tstp+1] = rhoN3tilde.block(6,6,10,10).trace();
      EXCN3_2_exp[tstp+1] = rhoN3tilde.block(0,0,6,6).trace();
      Eigen::ComplexEigenSolver<Eigen::MatrixXcd> eS3(rhoN3,false);
      eigN3 = eS3.eigenvalues();
      //
      // N=4 sector
      rhoN4.block(0,0,3,3)  =rho[tstp].M_[7] ;
      rhoN4.block(3,3,9,9)  =rho[tstp].M_[10];
      rhoN4.block(12,12,3,3)=rho[tstp].M_[13];
      rhoN4tilde=rotdagN4*rhoN4*rotN4;
      GSN4_exp[tstp+1]    = rhoN4tilde.block(6,6,9,9).trace();
      EXCN4_1_exp[tstp+1] = rhoN4tilde.block(1,1,5,5).trace();
      EXCN4_2_exp[tstp+1] = rhoN4tilde.block(0,0,1,1).trace();
      Eigen::ComplexEigenSolver<Eigen::MatrixXcd> eS4(rhoN4,false);
      eigN4 = eS4.eigenvalues();
      //
      // N=5 sector
      rhoN5.block(0,0,3,3)  =rho[tstp].M_[11];
      rhoN5.block(3,3,3,3)  =rho[tstp].M_[14];
      GSN5_exp[tstp+1] = rhoN5.trace();
      //
      // N=6 sector
      rhoN6.block(0,0,1,1)  =rho[tstp].M_[15];
    }

    //
    std::cout << " tstp, n, n1, n2, n3, d1, d2 , d3, "
              <<tstp<<", "
              <<n_exp[tstp+1]<<", " <<n1_exp[tstp+1]<<", "<<n2_exp[tstp+1]<<", "<<n3_exp[tstp+1]<<", "
              <<d1_exp[tstp+1]<<", "<<d2_exp[tstp+1]<<", "<<d3_exp[tstp+1]<<", "
              <<std::endl;

    if(tstp>5)
    {
      //
      std::string pathdir="./picdir_s";
      std::string sitedir=pathdir.append(std::to_string(site+1)).append("/");
      //
      std::string f1  = sitedir+"n.out"           ; const char * File1   = f1.c_str() ; FILE * outFile1 ; outFile1  = fopen(File1 ,"a");
      std::string f2  = sitedir+"d.out"           ; const char * File2   = f2.c_str() ; FILE * outFile2 ; outFile2  = fopen(File2 ,"a");
      std::string f3  = sitedir+"m.out"           ; const char * File3   = f3.c_str() ; FILE * outFile3 ; outFile3  = fopen(File3 ,"a");
      std::string f4  = sitedir+"nn_intra.out"    ; const char * File4   = f4.c_str() ; FILE * outFile4 ; outFile4  = fopen(File4 ,"a");
      std::string f5  = sitedir+"nn_inter.out"    ; const char * File5   = f5.c_str() ; FILE * outFile5 ; outFile5  = fopen(File5 ,"a");
      std::string f6  = sitedir+"nn_inter_eqS.out"; const char * File6   = f6.c_str() ; FILE * outFile6 ; outFile6  = fopen(File6 ,"a");
      std::string f7  = sitedir+"nnnn_ints.out"   ; const char * File7   = f7.c_str() ; FILE * outFile7 ; outFile7  = fopen(File7 ,"a");
      std::string f8  = sitedir+"Epot.out"        ; const char * File8   = f8.c_str() ; FILE * outFile8 ; outFile8  = fopen(File8 ,"a");
      std::string f9  = sitedir+"rho.out"         ; const char * File9   = f9.c_str() ; FILE * outFile9 ; outFile9  = fopen(File9 ,"a");
      //
      std::string f10 = sitedir+"LocStatesN0.out" ; const char * File10  = f10.c_str(); FILE * outFile10; outFile10 = fopen(File10,"a");
      std::string f11 = sitedir+"LocStatesN1.out" ; const char * File11  = f11.c_str(); FILE * outFile11; outFile11 = fopen(File11,"a");
      std::string f12 = sitedir+"LocStatesN2.out" ; const char * File12  = f12.c_str(); FILE * outFile12; outFile12 = fopen(File12,"a");
      std::string f13 = sitedir+"LocStatesN3.out" ; const char * File13  = f13.c_str(); FILE * outFile13; outFile13 = fopen(File13,"a");
      std::string f14 = sitedir+"LocStatesN4.out" ; const char * File14  = f14.c_str(); FILE * outFile14; outFile14 = fopen(File14,"a");
      std::string f15 = sitedir+"LocStatesN5.out" ; const char * File15  = f15.c_str(); FILE * outFile15; outFile15 = fopen(File15,"a");
      std::string f16 = sitedir+"LocStatesN6.out" ; const char * File16  = f16.c_str(); FILE * outFile16; outFile16 = fopen(File16,"a");
      //
      std::string f17 = sitedir+"SumRule.out"     ; const char * File17  = f17.c_str(); FILE * outFile17; outFile17 = fopen(File17,"a");
      //
      std::string f18 = sitedir+"EigenN2.out"     ; const char * File18  = f18.c_str(); FILE * outFile18; outFile18 = fopen(File18,"a");
      std::string f19 = sitedir+"EigenN3.out"     ; const char * File19  = f19.c_str(); FILE * outFile19; outFile19 = fopen(File19,"a");
      std::string f20 = sitedir+"EigenN4.out"     ; const char * File20  = f20.c_str(); FILE * outFile20; outFile20 = fopen(File20,"a");
      //
      fprintf (outFile1 , "%i\t%.20e\t%.20e\t%.20e\t%.20e\n",tstp,n1_exp[tstp+1],n2_exp[tstp+1],n3_exp[tstp+1],n_exp[tstp+1]);                                 fclose(outFile1);
      fprintf (outFile2 , "%i\t%.20e\t%.20e\t%.20e\t%.20e\n",tstp,d1_exp[tstp+1],d2_exp[tstp+1],d3_exp[tstp+1],(d1_exp[tstp+1]+d2_exp[tstp+1]+d3_exp[tstp+1]));fclose(outFile2);
      fprintf (outFile3 , "%i\t%.20e\t%.20e\t%.20e\t%.20e\n",tstp,m1_exp[tstp+1],m2_exp[tstp+1],m3_exp[tstp+1],m_exp[tstp+1]);                                 fclose(outFile3);
      fprintf (outFile4 , "%i\t%.20e\n",tstp,densdens_intraband_exp[tstp+1]);                                                                                  fclose(outFile4);
      fprintf (outFile5 , "%i\t%.20e\n",tstp,densdens_interband_exp[tstp+1]);                                                                                  fclose(outFile5);
      fprintf (outFile6 , "%i\t%.20e\n",tstp,densdens_interband_equalspin_exp[tstp+1]);                                                                        fclose(outFile6);
      fprintf (outFile7 , "%i\t%.20e\t%.20e\t%.20e\n",tstp,spin_flip_exp[tstp+1],pair_hop_exp[tstp+1],Eint_exp[tstp+1]);                                       fclose(outFile7);
      fprintf (outFile8 , "%i\t%.20e\n",tstp,Eint_exp[tstp+1]);                                                                                                fclose(outFile8);
      fprintf (outFile9 , "%i\t%.20e\t%.20e\n",tstp, expectation_value(tstp,rho,Q).real(), expectation_value(tstp,rho,Q).imag());                              fclose(outFile9);
      //
      fprintf (outFile10, "%i\t%.20e\t%.20e\n",tstp, rhoN0(0,0).real(), rhoN0(0,0).imag());                                                                    fclose(outFile10);
      fprintf (outFile11, "%i\t%.20e\t%.20e\n",tstp,GSN1_exp[tstp+1].real(),GSN1_exp[tstp+1].imag());                                                          fclose(outFile11);
      fprintf (outFile12, "%i\t%.20e\t%.20e\t%.20e\t%.20e\t%.20e\t%.20e\n",tstp,GSN2_exp[tstp+1].real(),GSN2_exp[tstp+1].imag(),
      EXCN2_1_exp[tstp+1].real(),EXCN2_1_exp[tstp+1].imag(),EXCN2_2_exp[tstp+1].real(),EXCN2_2_exp[tstp+1].imag());                                            fclose(outFile12);
      fprintf (outFile13, "%i\t%.20e\t%.20e\t%.20e\t%.20e\t%.20e\t%.20e\n",tstp,GSN3_exp[tstp+1].real(),GSN3_exp[tstp+1].imag(),
      EXCN3_1_exp[tstp+1].real(),EXCN3_1_exp[tstp+1].imag(),EXCN3_2_exp[tstp+1].real(),EXCN3_2_exp[tstp+1].imag());                                            fclose(outFile13);
      fprintf (outFile14, "%i\t%.20e\t%.20e\t%.20e\t%.20e\t%.20e\t%.20e\n",tstp,GSN4_exp[tstp+1].real(),GSN4_exp[tstp+1].imag(),
      EXCN4_1_exp[tstp+1].real(),EXCN4_1_exp[tstp+1].imag(),EXCN4_2_exp[tstp+1].real(),EXCN4_2_exp[tstp+1].imag());                                            fclose(outFile14);
      fprintf (outFile15, "%i\t%.20e\t%.20e\n",tstp,GSN5_exp[tstp+1].real(),GSN5_exp[tstp+1].imag());                                                          fclose(outFile15);
      fprintf (outFile16, "%i\t%.20e\t%.20e\n",tstp, rhoN6(0,0).real(), rhoN6(0,0).imag());                                                                    fclose(outFile16);
      //
      fprintf (outFile17, "%i\t%.20e\t%.20e\t%.20e\t%.20e\t%.20e\t%.20e\t%.20e\t%.20e\t%.20e\t%.20e\t%.20e\t%.20e\t%.20e\t%.20e\t%.20e\t%.20e\t%.20e\t%.20e\t%.20e\t%.20e\t%.20e\t%.20e\t%.20e\n",
                                               tstp,rhoN0(0,0).real(),rho[tstp].M_[0].trace().real()
                                                   ,rhoN1.trace().real(),GSN1_exp[tstp+1].real()
                                                   ,rhoN2.trace().real(),rhoN2tilde.trace().real(),GSN2_exp[tstp+1].real(),EXCN2_1_exp[tstp+1].real(),EXCN2_2_exp[tstp+1].real()
                                                   ,rhoN3.trace().real(),rhoN3tilde.trace().real(),GSN3_exp[tstp+1].real(),EXCN3_1_exp[tstp+1].real(),EXCN3_2_exp[tstp+1].real()
                                                   ,rhoN4.trace().real(),rhoN4tilde.trace().real(),GSN4_exp[tstp+1].real(),EXCN4_1_exp[tstp+1].real(),EXCN4_2_exp[tstp+1].real()
                                                   ,rhoN5.trace().real(),GSN5_exp[tstp+1].real()
                                                   ,rhoN6(0,0).real(),rho[tstp].M_[15].trace().real());                                                        fclose(outFile17);
      //
      fprintf (outFile18, "%i\t",tstp);
      for(int ie=0; ie < 15; ie++)fprintf (outFile18, "%.20e\t",eigN2(ie).real());
      fprintf (outFile18, "\n");fclose(outFile18);
      //
      fprintf (outFile19, "%i\t",tstp);
      for(int ie=0; ie < 20; ie++)fprintf (outFile19, "%.20e\t",eigN3(ie).real());
      fprintf (outFile19, "\n");fclose(outFile19);
      //
      fprintf (outFile20, "%i\t",tstp);
      for(int ie=0; ie < 15; ie++)fprintf (outFile20, "%.20e\t",eigN4(ie).real());
      fprintf (outFile20, "\n");fclose(outFile20);

   }

   if(tstp==10)
   {
      std::cout << "PRINTING OF RHO"  << std::endl;
      std::string pathdir="./picdir_s";
      std::string sitedir=pathdir.append(std::to_string(site+1)).append("/");
      std::ofstream file(sitedir+"MATRIX.out");
      //
      file << "##########   N=0  ##########" << '\n';
      file << "rho.M_ " << '\n';
      file << rho[tstp].M_[0].real().format(myFmt) << '\n'<< '\n';
      file << "joined rhoN " << '\n';
      file << rhoN0.real().format(myFmt) << '\n'<< '\n';
      //
      file << '\n'<< '\n'<< '\n';
      file << "##########   N=1  ##########" << '\n';
      file << "rho.M_ " << '\n';
      file << rho[tstp].M_[1].real().format(myFmt) << '\n'<< '\n';
      file << rho[tstp].M_[4].real().format(myFmt) << '\n'<< '\n';
      file << "joined rhoN " << '\n';
      file << rhoN1.real().format(myFmt) << '\n'<< '\n';
      //
      file << '\n'<< '\n'<< '\n';
      file << "##########   N=2  ##########" << '\n';
      file << "rho.M_ " << '\n';
      file << rho[tstp].M_[2].real().format(myFmt) << '\n'<< '\n';
      file << rho[tstp].M_[5].real().format(myFmt) << '\n'<< '\n';
      file << rho[tstp].M_[8].real().format(myFmt) << '\n'<< '\n';
      file << "joined rhoN " << '\n';
      file << rhoN2.real().format(myFmt) << '\n'<< '\n';
      file << "rotated rhoN " << '\n';
      file << rhoN2tilde.real().format(myFmt) << '\n'<< '\n';
      //
      file << '\n'<< '\n'<< '\n';
      file << "##########   N=3  ##########" << '\n';
      file << "rho.M_ " << '\n';
      file << rho[tstp].M_[3].real().format(myFmt) << '\n'<< '\n';
      file << rho[tstp].M_[6].real().format(myFmt) << '\n'<< '\n';
      file << rho[tstp].M_[9].real().format(myFmt) << '\n'<< '\n';
      file << rho[tstp].M_[12].real().format(myFmt) << '\n'<< '\n';
      file << "joined rhoN " << '\n';
      file << rhoN3.real().format(myFmt) << '\n'<< '\n';
      file << "rotated rhoN " << '\n';
      file << rhoN3tilde.real().format(myFmt) << '\n'<< '\n';
      //
      file << '\n'<< '\n'<< '\n';
      file << "##########   N=4  ##########" << '\n';
      file << "rho.M_ " << '\n';
      file << rho[tstp].M_[7].real().format(myFmt) << '\n'<< '\n';
      file << rho[tstp].M_[10].real().format(myFmt) << '\n'<< '\n';
      file << rho[tstp].M_[13].real().format(myFmt) << '\n'<< '\n';
      file << "joined rhoN " << '\n';
      file << rhoN4.real().format(myFmt) << '\n'<< '\n';
      file << "rotated rhoN " << '\n';
      file << rhoN4tilde.real().format(myFmt) << '\n'<< '\n';
      //
      file << '\n'<< '\n'<< '\n';
      file << "##########   N=5  ##########" << '\n';
      file << "rho.M_ " << '\n';
      file << rho[tstp].M_[11].real().format(myFmt) << '\n'<< '\n';
      file << rho[tstp].M_[14].real().format(myFmt) << '\n'<< '\n';
      file << "joined rhoN " << '\n';
      file << rhoN5.real().format(myFmt) << '\n'<< '\n';
      //
      file << '\n'<< '\n'<< '\n';
      file << "##########   N=6  ##########" << '\n';
      file << "rho.M_ " << '\n';
      file << rho[tstp].M_[15].real().format(myFmt) << '\n'<< '\n';
      file << "joined rhoN " << '\n';
      file << rhoN6.real().format(myFmt) << '\n'<< '\n'<< '\n'<< '\n'<< '\n'<< '\n'<< '\n'<< '\n';
      file.close();
   }
  }

  // ---------------------------------------------------------------------
  void store(hid_t group_id) {
    //
    store_double_attribute_to_hid(group_id, "mu", mu);
    store_real_data_to_hid(group_id, "U", U.data(), U.size());
    store_real_data_to_hid(group_id, "J", J.data(), J.size());
    //
    store_real_data_to_hid(group_id, "Q_exp"   , Q_exp.data(), Q_exp.size());
    store_real_data_to_hid(group_id, "Eint_exp", Eint_exp.data(), Eint_exp.size());
    //
    store_real_data_to_hid(group_id, "n_exp" , n_exp.data() , n_exp.size());
    store_real_data_to_hid(group_id, "n1_exp", n1_exp.data(), n1_exp.size());
    store_real_data_to_hid(group_id, "n2_exp", n2_exp.data(), n2_exp.size());
    store_real_data_to_hid(group_id, "n3_exp", n3_exp.data(), n3_exp.size());
    //
    store_real_data_to_hid(group_id, "m_exp" , m_exp.data() , m_exp.size());
    store_real_data_to_hid(group_id, "m1_exp", m1_exp.data(), m1_exp.size());
    store_real_data_to_hid(group_id, "m2_exp", m2_exp.data(), m2_exp.size());
    store_real_data_to_hid(group_id, "m3_exp", m3_exp.data(), m3_exp.size());
    //
    store_real_data_to_hid(group_id, "d1_exp", d1_exp.data(), d1_exp.size());
    store_real_data_to_hid(group_id, "d2_exp", d2_exp.data(), d2_exp.size());
    store_real_data_to_hid(group_id, "d3_exp", d3_exp.data(), d3_exp.size());
    //
    store_real_data_to_hid(group_id, "densdens_intraband_exp"           , densdens_intraband_exp.data(), densdens_intraband_exp.size());
    store_real_data_to_hid(group_id, "densdens_interband_exp"           , densdens_interband_exp.data(), densdens_interband_exp.size());
    store_real_data_to_hid(group_id, "densdens_interband_equalspin_exp" , densdens_interband_equalspin_exp.data(), densdens_interband_equalspin_exp.size());
    store_real_data_to_hid(group_id, "spin_flip_exp"                    , spin_flip_exp.data(), spin_flip_exp.size());
    store_real_data_to_hid(group_id, "pair_hop_exp"                     , pair_hop_exp.data(), pair_hop_exp.size());
    //
    store_cplx_data_to_hid(group_id, "GSN1_exp"      , GSN1_exp.data()      , GSN1_exp.size()      );
    store_cplx_data_to_hid(group_id, "GSN2_exp"      , GSN2_exp.data()      , GSN2_exp.size()      );
    store_cplx_data_to_hid(group_id, "EXCN2_1_exp"   , EXCN2_1_exp.data()   , EXCN2_1_exp.size()   );
    store_cplx_data_to_hid(group_id, "EXCN2_2_exp"   , EXCN2_2_exp.data()   , EXCN2_2_exp.size()   );
    store_cplx_data_to_hid(group_id, "GSN3_exp"      , GSN3_exp.data()      , GSN3_exp.size()      );
    store_cplx_data_to_hid(group_id, "EXCN3_1_exp"   , EXCN3_1_exp.data()   , EXCN3_1_exp.size()   );
    store_cplx_data_to_hid(group_id, "EXCN3_2_exp"   , EXCN3_2_exp.data()   , EXCN3_2_exp.size()   );
    store_cplx_data_to_hid(group_id, "GSN4_exp"      , GSN4_exp.data()      , GSN4_exp.size()      );
    store_cplx_data_to_hid(group_id, "EXCN4_1_exp"   , EXCN4_1_exp.data()   , EXCN4_1_exp.size()   );
    store_cplx_data_to_hid(group_id, "EXCN4_2_exp"   , EXCN4_2_exp.data()   , EXCN4_2_exp.size()   );
    store_cplx_data_to_hid(group_id, "GSN5_exp"      , GSN5_exp.data()      , GSN5_exp.size()      );
    //

    // -- Sub group operators:
    hid_t sub_group_id = create_group(group_id, "op");
    store_operator_type(sub_group_id, c1uc, std::string("c1uc"));
    store_operator_type(sub_group_id, c1ua, std::string("c1ua"));
    store_operator_type(sub_group_id, c1dc, std::string("c1dc"));
    store_operator_type(sub_group_id, c1da, std::string("c1da"));
    store_operator_type(sub_group_id, c2uc, std::string("c2uc"));
    store_operator_type(sub_group_id, c2ua, std::string("c2ua"));
    store_operator_type(sub_group_id, c2dc, std::string("c2dc"));
    store_operator_type(sub_group_id, c2da, std::string("c2da"));
    store_operator_type(sub_group_id, c3uc, std::string("c3uc"));
    store_operator_type(sub_group_id, c3ua, std::string("c3ua"));
    store_operator_type(sub_group_id, c3dc, std::string("c3dc"));
    store_operator_type(sub_group_id, c3da, std::string("c3da"));
    store_operator_type(sub_group_id, n1u, std::string("n1u"));
    store_operator_type(sub_group_id, n1d, std::string("n1d"));
    store_operator_type(sub_group_id, n2u, std::string("n2u"));
    store_operator_type(sub_group_id, n2d, std::string("n2d"));
    store_operator_type(sub_group_id, n3u, std::string("n3u"));
    store_operator_type(sub_group_id, n3d, std::string("n3d"));
    store_operator_type(sub_group_id, n , std::string("n"));
    store_operator_type(sub_group_id, n1, std::string("n1"));
    store_operator_type(sub_group_id, n2, std::string("n2"));
    store_operator_type(sub_group_id, n3, std::string("n3"));
    store_operator_type(sub_group_id, m , std::string("m"));
    store_operator_type(sub_group_id, m1, std::string("m1"));
    store_operator_type(sub_group_id, m2, std::string("m2"));
    store_operator_type(sub_group_id, m3, std::string("m3"));
    store_operator_type(sub_group_id, d1, std::string("d1"));
    store_operator_type(sub_group_id, d2, std::string("d2"));
    store_operator_type(sub_group_id, d3, std::string("d3"));
    store_operator_type(sub_group_id, Q, std::string("Q"));

    store_operator_type(sub_group_id, densdens_intraband, std::string("densdens_intraband"));
    store_operator_type(sub_group_id, densdens_interband, std::string("densdens_interband"));
    store_operator_type(sub_group_id, densdens_interband_equalspin, std::string("densdens_interband_equalspin"));

    store_operator_type(sub_group_id, spin_flip, std::string("spin_flip"));
    store_operator_type(sub_group_id, pair_hop, std::string("pair_hop"));

    close_group(sub_group_id);

  }

  int nt;
  HILB hil_;
  //
  int site;
  double mu;
  std::vector<double> U, J, Bz;
  std::vector<double> E1, E2, E3;
  std::vector<double> Q_exp, Eint_exp;
  std::vector<double> n1_exp, n2_exp, n3_exp, n_exp;
  std::vector<double> m1_exp, m2_exp, m3_exp, m_exp;
  std::vector<double> d1_exp, d2_exp, d3_exp;
  //
  std::vector<double> densdens_intraband_exp;
  std::vector<double> densdens_interband_exp;
  std::vector<double> densdens_interband_equalspin_exp;
  std::vector<double> spin_flip_exp, pair_hop_exp;
  //
  std::vector<std::complex<double>> GSN1_exp;
  std::vector<std::complex<double>> GSN2_exp, EXCN2_1_exp, EXCN2_2_exp;
  std::vector<std::complex<double>> GSN3_exp, EXCN3_1_exp, EXCN3_2_exp;
  std::vector<std::complex<double>> GSN4_exp, EXCN4_1_exp, EXCN4_2_exp;
  std::vector<std::complex<double>> GSN5_exp;
  //
  operator_type c1ua, c1uc, c1da, c1dc,	c2ua, c2uc, c2da, c2dc, c3ua, c3uc, c3da, c3dc;
  //
  operator_type n1u, n1d, n2u, n2d, n3u, n3d;
  operator_type n1, n2, n3, n, Q;
  operator_type m1, m2, m3, m;
  operator_type d1, d2, d3;
  //
  operator_type densdens_intraband;
  operator_type densdens_interband;
  operator_type densdens_interband_equalspin;
  operator_type spin_flip, pair_hop;
  //
  Eigen::MatrixXcd rotN2,rotN3,rotN4,rotdagN2,rotdagN3,rotdagN4;

};

// -----------------------------------------------------------------------
} // end namespace hamiltonians
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
} // end namespace ppsc
// -----------------------------------------------------------------------

#endif // _PPSC_HAMILTONIANS_THREE_BAND_HUBBARD_SOLAR_HPP







/*
if(tstp==2)
{
Eigen::IOFormat myFmt(3, 0, ", ", ";\n", "", "", "[", "]");
for(int s=0; s < rho[tstp].ns_; s++)
  {
    std::cout << "begin de la mona de tu mare vaca"  << std::endl;
    std::cout<< "sector: "<< s<<"  sector dim: "<< rho[tstp].ssdim_[s] << std::endl;
    std::cout<< "rho" <<std::endl;
    std::cout<< rho[tstp].M_[s].format(myFmt)<< std::endl;
    std::cout << "end de la mona de tu mare vaca"  << std::endl;
  }
  Eigen::MatrixXcd rhoN1(6,6)   ;   rhoN1.setZero(6,6) ; rhoN1.block(0,0,3,3)=rho[tstp].M_[1] ; rhoN1.block(3,3,3,3)=rho[tstp].M_[4] ;
  Eigen::MatrixXcd A(15,15), rhoN2(15,15) ; rhoN2.setZero(15,15) ; rhoN2.block(0,0,3,3)=rho[tstp].M_[2] ; rhoN2.block(3,3,9,9)=rho[tstp].M_[5] ; rhoN2.block(12,12,3,3)=rho[tstp].M_[8] ;
  Eigen::MatrixXcd B(20,20), rhoN3(20,20) ; rhoN3.setZero(20,20) ; rhoN3.block(0,0,1,1)=rho[tstp].M_[3] ; rhoN3.block(1,1,9,9)=rho[tstp].M_[6] ; rhoN3.block(10,10,9,9)=rho[tstp].M_[9] ; rhoN3.block(19,19,1,1)=rho[tstp].M_[12] ;

  std::cout << " "  << std::endl;std::cout << " "  << std::endl;std::cout << " "  << std::endl;
  std::cout << "////////////////////////////////////////////////////////"  << std::endl;
  std::cout << " "  << std::endl;
  std::cout << "rho N=1"  << std::endl;
  std::cout<<rhoN1.format(myFmt)<<std::endl;
  std::cout << "////////////////////////////////////////////////////////"  << std::endl;

  std::cout << " "  << std::endl;std::cout << " "  << std::endl;std::cout << " "  << std::endl;
  std::cout << "////////////////////////////////////////////////////////"  << std::endl;
  std::cout << " "  << std::endl;
  std::cout << "rho N=2"  << std::endl;
  std::cout<<rhoN2.format(myFmt)<<std::endl;
  std::cout << " "  << std::endl;
  std::cout << "rotrho N=2"  << std::endl;
  A=rotdagN2*rhoN2*rotN2;
  std::cout<<A.format(myFmt)<<std::endl;
  std::cout<<"        "<<std::endl;
  std::cout<<A.real().format(myFmt)<<std::endl;

  std::cout << "////////////////////////////////////////////////////////"  << std::endl;


  std::cout << " "  << std::endl;std::cout << " "  << std::endl;std::cout << " "  << std::endl;
  std::cout << "////////////////////////////////////////////////////////"  << std::endl;
  std::cout << " "  << std::endl;
  std::cout << "rho N=3"  << std::endl;
  std::cout<<rhoN3.format(myFmt)<<std::endl;
  std::cout << " "  << std::endl;
  std::cout << "rotrho N=3"  << std::endl;
  B=rotdagN3*rhoN3*rotN3;
  std::cout<<B.format(myFmt)<<std::endl;
  std::cout << "////////////////////////////////////////////////////////"  << std::endl;
}

*/


/*

operator_type GSN1;
operator_type GSN2_Szp1, GSN2_Sz0, GSN2_Szm1, EXCN2_1, EXCN2_2;
operator_type GSN3_Szp32, GSN3_Szp12, GSN3_Szm12, GSN3_Szm32, EXCN3_1, EXCN3_2;

double oversqrttwo = 1/sqrt(2);
double oversqrtthree = 1/sqrt(3);
// [N=1]
// Gs (x6)
GSN1 = (n1u) + (n2u) + (n3u) + (n1d) + (n2d) + (n3d);
// [N=2]
// Gs (x9)
GSN2_Szp1 = (n2u*n3u) + (n1u*n3u) + (n1u*n2u);
GSN2_Sz0  = (n2u*n3d+n2d*n3u)*oversqrttwo + (n1u*n3d+n1d*n3u)*oversqrttwo + (n1u*n2d+n1d*n2u)*oversqrttwo;
GSN2_Szm1 = (n2d*n3d) + (n1d*n3d) + (n1d*n2d);
// Exctd-1 (x5)
EXCN2_1 = (n1u*n1d-n3u*n3d)*oversqrttwo + (n2u*n2d-n3u*n3d)*oversqrttwo
        + (n2u*n3d-n2d*n3u)*oversqrttwo + (n1d*n2u-n1u*n2d)*oversqrttwo
        + (n1u*n3d-n1d*n3u)*oversqrttwo;
// Exctd-2 (x1)
EXCN2_2 = (n1u*n1d+n2u*n2d+n3u*n3d)*oversqrtthree;
// [N=3]
// Gs (x9)
GSN3_Szp32 = (n1u*n2u*n3u);
GSN3_Szp12 = (n1d*n2u*n3u+n1u*n2d*n3u+n1u*n2u*n3d)*oversqrtthree;
GSN3_Szm12 = (n1d*n2d*n3u+n1d*n2u*n3d+n1u*n2d*n3d)*oversqrtthree;
GSN3_Szm32 = (n1d*n2d*n3d);
// Exctd-1 (x10)
EXCN3_1 = (n1u*n1d*n2u-n2u*n3u*n3d)*oversqrttwo + (n1u*n2u*n2d-n1u*n3u*n3d)*oversqrttwo
        + (n1u*n1d*n2d-n2d*n3u*n3d)*oversqrttwo + (n1u*n2u*n3d-n1d*n2u*n3u)*oversqrttwo
        + (n1u*n1d*n3u-n2u*n2d*n3u)*oversqrttwo + (n1u*n2d*n3u-n1d*n2u*n3u)*oversqrttwo
        + (n1u*n1d*n3d-n2u*n2d*n3d)*oversqrttwo + (n1u*n2d*n3d-n1d*n2d*n3u)*oversqrttwo
        + (n1d*n2u*n2d-n1d*n3u*n3d)*oversqrttwo + (n1d*n2u*n3d-n1d*n2d*n3u)*oversqrttwo;
// Exctd-2 (x6)
EXCN3_2 = (n2u*n3u*n3d+n1u*n1d*n2u)*oversqrttwo + (n2d*n3u*n3d+n1u*n1d*n2d)*oversqrttwo
        + (n2u*n2d*n3u+n1u*n1d*n3u)*oversqrttwo + (n2u*n2d*n3d+n1u*n1d*n3d)*oversqrttwo
        + (n1u*n3u*n3d+n1u*n2u*n2d)*oversqrttwo + (n1d*n3u*n3d+n1d*n2u*n2d)*oversqrttwo;
*/

// -----------------------------------------------------------------------------
// This file was autogenerated by symforce from template:
//     function/FUNCTION.h.jinja
// Do NOT modify by hand.
// -----------------------------------------------------------------------------

#pragma once

#include <Eigen/Dense>

#include <sym/rot3.h>

namespace sym {

/**
 * This function was autogenerated from a symbolic function. Do not modify by hand.
 *
 * Symbolic function: IK_residual_func_cost1_Nl7
 *
 * Args:
 *     fh1: Scalar
 *     fv1: Scalar
 *     DeltaRot: Rot3
 *     p_init0: Scalar
 *     p_init1: Scalar
 *     p_init2: Scalar
 *     rot_init_x: Scalar
 *     rot_init_y: Scalar
 *     rot_init_z: Scalar
 *     rot_init_w: Scalar
 *     epsilon: Scalar
 *
 * Outputs:
 *     res: Matrix41
 */
template <typename Scalar>
Eigen::Matrix<Scalar, 4, 1> IkResidualFuncCost1Nl7(const Scalar fh1, const Scalar fv1,
                                                   const sym::Rot3<Scalar>& DeltaRot,
                                                   const Scalar p_init0, const Scalar p_init1,
                                                   const Scalar p_init2, const Scalar rot_init_x,
                                                   const Scalar rot_init_y, const Scalar rot_init_z,
                                                   const Scalar rot_init_w, const Scalar epsilon) {
  // Total ops: 507

  // Unused inputs
  (void)epsilon;

  // Input arrays
  const Eigen::Matrix<Scalar, 4, 1>& _DeltaRot = DeltaRot.Data();

  // Intermediate terms (156)
  const Scalar _tmp0 = Scalar(1.0) / (fh1);
  const Scalar _tmp1 = _DeltaRot[0] * rot_init_z + _DeltaRot[1] * rot_init_w -
                       _DeltaRot[2] * rot_init_x + _DeltaRot[3] * rot_init_y;
  const Scalar _tmp2 = -2 * std::pow(_tmp1, Scalar(2));
  const Scalar _tmp3 = -_DeltaRot[0] * rot_init_y + _DeltaRot[1] * rot_init_x +
                       _DeltaRot[2] * rot_init_w + _DeltaRot[3] * rot_init_z;
  const Scalar _tmp4 = 1 - 2 * std::pow(_tmp3, Scalar(2));
  const Scalar _tmp5 = Scalar(0.20999999999999999) * _tmp2 + Scalar(0.20999999999999999) * _tmp4;
  const Scalar _tmp6 = _DeltaRot[0] * rot_init_w - _DeltaRot[1] * rot_init_z +
                       _DeltaRot[2] * rot_init_y + _DeltaRot[3] * rot_init_x;
  const Scalar _tmp7 = 2 * _tmp3 * _tmp6;
  const Scalar _tmp8 = -2 * _DeltaRot[0] * rot_init_x - 2 * _DeltaRot[1] * rot_init_y -
                       2 * _DeltaRot[2] * rot_init_z + 2 * _DeltaRot[3] * rot_init_w;
  const Scalar _tmp9 = _tmp1 * _tmp8;
  const Scalar _tmp10 = _tmp7 + _tmp9;
  const Scalar _tmp11 = -Scalar(0.010999999999999999) * _tmp10;
  const Scalar _tmp12 = 2 * _tmp1;
  const Scalar _tmp13 = _tmp12 * _tmp6;
  const Scalar _tmp14 = _tmp3 * _tmp8;
  const Scalar _tmp15 = Scalar(0.20999999999999999) * _tmp13 - Scalar(0.20999999999999999) * _tmp14;
  const Scalar _tmp16 = _tmp11 - _tmp15;
  const Scalar _tmp17 = _tmp16 + _tmp5;
  const Scalar _tmp18 = _tmp17 + p_init0;
  const Scalar _tmp19 = Scalar(0.20999999999999999) * _tmp13 + Scalar(0.20999999999999999) * _tmp14;
  const Scalar _tmp20 = _tmp12 * _tmp3;
  const Scalar _tmp21 = _tmp6 * _tmp8;
  const Scalar _tmp22 = _tmp20 - _tmp21;
  const Scalar _tmp23 = -Scalar(0.010999999999999999) * _tmp22;
  const Scalar _tmp24 = -2 * std::pow(_tmp6, Scalar(2));
  const Scalar _tmp25 = Scalar(0.20999999999999999) * _tmp24 + Scalar(0.20999999999999999) * _tmp4;
  const Scalar _tmp26 = _tmp23 - _tmp25;
  const Scalar _tmp27 = _tmp19 + _tmp26;
  const Scalar _tmp28 = _tmp27 + p_init1;
  const Scalar _tmp29 = std::asinh(_tmp0 * fv1);
  const Scalar _tmp30 = Scalar(9.6622558468725703) * fh1;
  const Scalar _tmp31 = Scalar(0.20999999999999999) * _tmp20 + Scalar(0.20999999999999999) * _tmp21;
  const Scalar _tmp32 = -_tmp31;
  const Scalar _tmp33 = -Scalar(0.010999999999999999) * _tmp2 -
                        Scalar(0.010999999999999999) * _tmp24 + Scalar(-0.010999999999999999);
  const Scalar _tmp34 = Scalar(0.20999999999999999) * _tmp7 - Scalar(0.20999999999999999) * _tmp9;
  const Scalar _tmp35 = _tmp33 + _tmp34;
  const Scalar _tmp36 = _tmp32 + _tmp35;
  const Scalar _tmp37 = -_tmp19;
  const Scalar _tmp38 = _tmp23 + _tmp25;
  const Scalar _tmp39 = _tmp37 + _tmp38;
  const Scalar _tmp40 = _tmp39 + p_init1;
  const Scalar _tmp41 = _tmp40 + Scalar(-4.8333311099999996);
  const Scalar _tmp42 = -_tmp5;
  const Scalar _tmp43 = _tmp11 + _tmp15;
  const Scalar _tmp44 = _tmp42 + _tmp43;
  const Scalar _tmp45 = _tmp44 + p_init0;
  const Scalar _tmp46 = _tmp45 + Scalar(1.79662371);
  const Scalar _tmp47 = std::pow(Scalar(std::pow(_tmp41, Scalar(2)) + std::pow(_tmp46, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp48 = _tmp46 * _tmp47;
  const Scalar _tmp49 = _tmp43 + _tmp5;
  const Scalar _tmp50 = _tmp19 + _tmp38;
  const Scalar _tmp51 = _tmp50 + p_init1;
  const Scalar _tmp52 = _tmp51 + Scalar(-4.7752063900000001);
  const Scalar _tmp53 = _tmp49 + p_init0;
  const Scalar _tmp54 = _tmp53 + Scalar(-2.71799795);
  const Scalar _tmp55 = std::pow(Scalar(std::pow(_tmp52, Scalar(2)) + std::pow(_tmp54, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp56 = _tmp52 * _tmp55;
  const Scalar _tmp57 = _tmp54 * _tmp55;
  const Scalar _tmp58 = _tmp26 + _tmp37;
  const Scalar _tmp59 = _tmp58 + p_init1;
  const Scalar _tmp60 = _tmp59 + Scalar(8.3196563700000006);
  const Scalar _tmp61 = _tmp16 + _tmp42;
  const Scalar _tmp62 = _tmp61 + p_init0;
  const Scalar _tmp63 = _tmp62 + Scalar(1.9874742000000001);
  const Scalar _tmp64 =
      std::sqrt(Scalar(std::pow(_tmp60, Scalar(2)) + std::pow(_tmp63, Scalar(2))));
  const Scalar _tmp65 = Scalar(1.0) / (_tmp64);
  const Scalar _tmp66 = Scalar(1.0) / (_tmp63);
  const Scalar _tmp67 = _tmp64 * _tmp66;
  const Scalar _tmp68 = _tmp67 * (-_tmp58 * _tmp63 * _tmp65 + _tmp60 * _tmp61 * _tmp65);
  const Scalar _tmp69 = _tmp41 * _tmp47;
  const Scalar _tmp70 = _tmp39 * _tmp48 - _tmp44 * _tmp69 + _tmp48 * _tmp68;
  const Scalar _tmp71 = _tmp60 * _tmp66;
  const Scalar _tmp72 = Scalar(1.0) / (_tmp48 * _tmp71 - _tmp69);
  const Scalar _tmp73 = -_tmp56 + _tmp57 * _tmp71;
  const Scalar _tmp74 = _tmp72 * _tmp73;
  const Scalar _tmp75 = -_tmp49 * _tmp56 + _tmp50 * _tmp57 + _tmp57 * _tmp68 - _tmp70 * _tmp74;
  const Scalar _tmp76 = Scalar(1.0) / (_tmp75);
  const Scalar _tmp77 = Scalar(1.0) * _tmp72;
  const Scalar _tmp78 = _tmp73 * _tmp76 * _tmp77;
  const Scalar _tmp79 = Scalar(1.0) * _tmp76;
  const Scalar _tmp80 = _tmp18 + Scalar(-2.5202214700000001);
  const Scalar _tmp81 = _tmp28 + Scalar(8.3888750099999996);
  const Scalar _tmp82 = std::pow(Scalar(std::pow(_tmp80, Scalar(2)) + std::pow(_tmp81, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp83 = _tmp81 * _tmp82;
  const Scalar _tmp84 = _tmp80 * _tmp82;
  const Scalar _tmp85 = fh1 * (_tmp17 * _tmp83 - _tmp27 * _tmp84);
  const Scalar _tmp86 = Scalar(1.0) * _tmp61;
  const Scalar _tmp87 = -_tmp44 + _tmp86;
  const Scalar _tmp88 = _tmp33 - _tmp34;
  const Scalar _tmp89 = _tmp31 + _tmp88;
  const Scalar _tmp90 = _tmp32 + _tmp88;
  const Scalar _tmp91 = _tmp71 * _tmp90;
  const Scalar _tmp92 = -_tmp48 * _tmp91 + _tmp69 * _tmp89;
  const Scalar _tmp93 = Scalar(1.0) * _tmp58;
  const Scalar _tmp94 = -_tmp93;
  const Scalar _tmp95 = Scalar(1.0) / (_tmp39 + _tmp94);
  const Scalar _tmp96 = Scalar(1.0) * _tmp95;
  const Scalar _tmp97 = -_tmp48 * _tmp89 + _tmp48 * _tmp90;
  const Scalar _tmp98 = _tmp72 * _tmp87 * _tmp92 * _tmp96 - _tmp77 * _tmp97;
  const Scalar _tmp99 = _tmp31 + _tmp35;
  const Scalar _tmp100 = _tmp56 * _tmp99 - _tmp57 * _tmp91 - _tmp74 * _tmp92;
  const Scalar _tmp101 = _tmp87 * _tmp95;
  const Scalar _tmp102 = -_tmp100 * _tmp101 + _tmp57 * _tmp90 - _tmp57 * _tmp99 - _tmp74 * _tmp97;
  const Scalar _tmp103 = Scalar(1.0) / (_tmp102);
  const Scalar _tmp104 = _tmp103 * _tmp75;
  const Scalar _tmp105 = _tmp102 * _tmp76;
  const Scalar _tmp106 = _tmp105 * (-_tmp104 * _tmp98 - _tmp70 * _tmp77);
  const Scalar _tmp107 = _tmp106 + _tmp98;
  const Scalar _tmp108 = _tmp103 * _tmp73;
  const Scalar _tmp109 = -_tmp107 * _tmp108 + Scalar(1.0);
  const Scalar _tmp110 = _tmp48 * _tmp72;
  const Scalar _tmp111 = _tmp103 * _tmp57;
  const Scalar _tmp112 = _tmp83 * fh1;
  const Scalar _tmp113 = Scalar(43.164000000000001) - fv1;
  const Scalar _tmp114 = _tmp101 * _tmp93 + _tmp86;
  const Scalar _tmp115 = 0;
  const Scalar _tmp116 = _tmp108 * _tmp115;
  const Scalar _tmp117 = _tmp103 * _tmp115;
  const Scalar _tmp118 = _tmp71 * _tmp72;
  const Scalar _tmp119 = _tmp118 * _tmp92 + _tmp91;
  const Scalar _tmp120 = -_tmp101 * _tmp119 + _tmp118 * _tmp97 - _tmp90;
  const Scalar _tmp121 = _tmp105 * (-_tmp104 * _tmp120 + _tmp118 * _tmp70 - _tmp68);
  const Scalar _tmp122 = _tmp120 + _tmp121;
  const Scalar _tmp123 = -_tmp108 * _tmp122 - _tmp71;
  const Scalar _tmp124 = _tmp84 * fh1;
  const Scalar _tmp125 = -_tmp112 * _tmp67 * (_tmp107 * _tmp111 + _tmp109 * _tmp110) -
                         _tmp113 * _tmp67 * (-_tmp110 * _tmp116 + _tmp117 * _tmp57) -
                         _tmp124 * _tmp67 * (_tmp110 * _tmp123 + _tmp111 * _tmp122 + Scalar(1.0)) -
                         _tmp67 * _tmp85 * (-_tmp48 * _tmp78 + _tmp57 * _tmp79);
  const Scalar _tmp126 = Scalar(1.0) / (_tmp125);
  const Scalar _tmp127 = _tmp50 + _tmp94;
  const Scalar _tmp128 = _tmp101 * _tmp127;
  const Scalar _tmp129 = Scalar(1.0) / (-_tmp128 - _tmp49 + _tmp86);
  const Scalar _tmp130 = Scalar(1.0) * _tmp129;
  const Scalar _tmp131 = _tmp100 * _tmp103;
  const Scalar _tmp132 = _tmp127 * _tmp129;
  const Scalar _tmp133 = _tmp106 * _tmp132 - _tmp107 * _tmp131 - _tmp77 * _tmp92;
  const Scalar _tmp134 = _tmp128 * _tmp130 + Scalar(1.0);
  const Scalar _tmp135 = _tmp101 * _tmp130;
  const Scalar _tmp136 = _tmp36 * fh1;
  const Scalar _tmp137 = -_tmp136 * _tmp83 - Scalar(3.29616) * _tmp22 - _tmp27 * fv1;
  const Scalar _tmp138 = _tmp127 * _tmp95;
  const Scalar _tmp139 = Scalar(3.29616) * _tmp10 + _tmp136 * _tmp84 + _tmp17 * fv1;
  const Scalar _tmp140 = _tmp119 + _tmp121 * _tmp132 - _tmp122 * _tmp131;
  const Scalar _tmp141 = _tmp114 * _tmp129;
  const Scalar _tmp142 = _tmp95 * (-_tmp115 * _tmp131 - _tmp127 * _tmp141 + _tmp94);
  const Scalar _tmp143 = _tmp105 * _tmp130;
  const Scalar _tmp144 = -_tmp100 * _tmp79 + _tmp127 * _tmp143;
  const Scalar _tmp145 =
      std::asinh(_tmp126 * (Scalar(1.0) * _tmp112 * (_tmp106 * _tmp130 - _tmp133 * _tmp96) +
                            Scalar(1.0) * _tmp113 *
                                (-_tmp114 * _tmp130 - Scalar(1.0) * _tmp142 + Scalar(1.0)) +
                            Scalar(1.0) * _tmp124 * (_tmp121 * _tmp130 - _tmp140 * _tmp96) +
                            Scalar(1.0) * _tmp137 * (-_tmp134 * _tmp96 + _tmp135) +
                            Scalar(1.0) * _tmp139 * (_tmp130 * _tmp138 - _tmp130) +
                            Scalar(1.0) * _tmp85 * (_tmp143 - _tmp144 * _tmp96)));
  const Scalar _tmp146 = Scalar(9.6622558468725703) * _tmp125;
  const Scalar _tmp147 = _tmp130 * _tmp139;
  const Scalar _tmp148 = _tmp109 * _tmp112 * _tmp72 - _tmp113 * _tmp116 * _tmp72 +
                         _tmp123 * _tmp124 * _tmp72 - _tmp78 * _tmp85;
  const Scalar _tmp149 = Scalar(1.0) / (_tmp148);
  const Scalar _tmp150 = std::asinh(
      _tmp149 * (_tmp112 * _tmp133 * _tmp95 + _tmp113 * _tmp142 + _tmp124 * _tmp140 * _tmp95 +
                 _tmp134 * _tmp137 * _tmp95 - _tmp138 * _tmp147 + _tmp144 * _tmp85 * _tmp95));
  const Scalar _tmp151 = Scalar(9.6622558468725703) * _tmp148;
  const Scalar _tmp152 = _tmp103 * _tmp107 * _tmp112 + _tmp103 * _tmp122 * _tmp124 +
                         _tmp113 * _tmp117 + _tmp79 * _tmp85;
  const Scalar _tmp153 = Scalar(1.0) / (_tmp152);
  const Scalar _tmp154 = std::asinh(_tmp153 * (-_tmp106 * _tmp112 * _tmp129 + _tmp113 * _tmp141 -
                                               _tmp121 * _tmp124 * _tmp129 - _tmp135 * _tmp137 -
                                               _tmp143 * _tmp85 + _tmp147));
  const Scalar _tmp155 = Scalar(9.6622558468725703) * _tmp152;

  // Output terms (1)
  Eigen::Matrix<Scalar, 4, 1> _res;

  _res(0, 0) =
      -_tmp30 *
          (Scalar(0.87653584775870996) * _tmp0 + std::cosh(Scalar(1.0) * _tmp29) -
           std::cosh(
               Scalar(0.1034955) * _tmp0 *
               (-_tmp29 * _tmp30 -
                Scalar(8.3888750099999996) *
                    std::sqrt(Scalar(
                        Scalar(0.090254729040973036) *
                            std::pow(Scalar(1 - Scalar(0.39679052492160538) * _tmp18), Scalar(2)) +
                        std::pow(Scalar(-Scalar(0.11920549523123722) * _tmp28 - 1),
                                 Scalar(2))))))) +
      _tmp36 + p_init2;
  _res(1, 0) =
      -_tmp146 *
          (Scalar(0.87679799772039002) * _tmp126 + std::cosh(Scalar(1.0) * _tmp145) -
           std::cosh(
               Scalar(0.1034955) * _tmp126 *
               (-_tmp145 * _tmp146 -
                Scalar(8.3196563700000006) *
                    std::sqrt(Scalar(
                        std::pow(Scalar(-Scalar(0.12019727204189803) * _tmp59 - 1), Scalar(2)) +
                        Scalar(0.057067943376852184) *
                            std::pow(Scalar(-Scalar(0.50315118556004401) * _tmp62 - 1),
                                     Scalar(2))))))) +
      _tmp90 + p_init2;
  _res(2, 0) =
      -_tmp151 *
          (Scalar(0.86625939559540499) * _tmp149 + std::cosh(Scalar(1.0) * _tmp150) -
           std::cosh(
               Scalar(0.1034955) * _tmp149 *
               (-_tmp150 * _tmp151 -
                Scalar(4.8333311099999996) *
                    std::sqrt(Scalar(
                        std::pow(Scalar(1 - Scalar(0.20689664689659551) * _tmp40), Scalar(2)) +
                        Scalar(0.13817235445745474) *
                            std::pow(Scalar(-Scalar(0.55659957866191134) * _tmp45 - 1),
                                     Scalar(2))))))) +
      _tmp89 + p_init2;
  _res(3, 0) =
      -_tmp155 *
          (Scalar(0.86565325453551001) * _tmp153 + std::cosh(Scalar(1.0) * _tmp154) -
           std::cosh(
               Scalar(0.1034955) * _tmp153 *
               (-_tmp154 * _tmp155 -
                Scalar(4.7752063900000001) *
                    std::sqrt(Scalar(
                        std::pow(Scalar(1 - Scalar(0.20941503221602112) * _tmp51), Scalar(2)) +
                        Scalar(0.32397683292140877) *
                            std::pow(Scalar(1 - Scalar(0.36791786395571047) * _tmp53),
                                     Scalar(2))))))) +
      _tmp99 + p_init2;

  return _res;
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym

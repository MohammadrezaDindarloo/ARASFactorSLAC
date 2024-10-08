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
 * Symbolic function: IK_residual_func_cost1_Nl11
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
Eigen::Matrix<Scalar, 4, 1> IkResidualFuncCost1Nl11(
    const Scalar fh1, const Scalar fv1, const sym::Rot3<Scalar>& DeltaRot, const Scalar p_init0,
    const Scalar p_init1, const Scalar p_init2, const Scalar rot_init_x, const Scalar rot_init_y,
    const Scalar rot_init_z, const Scalar rot_init_w, const Scalar epsilon) {
  // Total ops: 499

  // Unused inputs
  (void)epsilon;

  // Input arrays
  const Eigen::Matrix<Scalar, 4, 1>& _DeltaRot = DeltaRot.Data();

  // Intermediate terms (153)
  const Scalar _tmp0 = Scalar(1.0) / (fh1);
  const Scalar _tmp1 = -_DeltaRot[0] * rot_init_y + _DeltaRot[1] * rot_init_x +
                       _DeltaRot[2] * rot_init_w + _DeltaRot[3] * rot_init_z;
  const Scalar _tmp2 = -2 * std::pow(_tmp1, Scalar(2));
  const Scalar _tmp3 = _DeltaRot[0] * rot_init_z + _DeltaRot[1] * rot_init_w -
                       _DeltaRot[2] * rot_init_x + _DeltaRot[3] * rot_init_y;
  const Scalar _tmp4 = 1 - 2 * std::pow(_tmp3, Scalar(2));
  const Scalar _tmp5 = Scalar(0.20999999999999999) * _tmp2 + Scalar(0.20999999999999999) * _tmp4;
  const Scalar _tmp6 = _DeltaRot[0] * rot_init_w - _DeltaRot[1] * rot_init_z +
                       _DeltaRot[2] * rot_init_y + _DeltaRot[3] * rot_init_x;
  const Scalar _tmp7 = 2 * _tmp1;
  const Scalar _tmp8 = _tmp6 * _tmp7;
  const Scalar _tmp9 = -2 * _DeltaRot[0] * rot_init_x - 2 * _DeltaRot[1] * rot_init_y -
                       2 * _DeltaRot[2] * rot_init_z + 2 * _DeltaRot[3] * rot_init_w;
  const Scalar _tmp10 = _tmp3 * _tmp9;
  const Scalar _tmp11 = _tmp10 + _tmp8;
  const Scalar _tmp12 = -Scalar(0.010999999999999999) * _tmp11;
  const Scalar _tmp13 = 2 * _tmp3 * _tmp6;
  const Scalar _tmp14 = _tmp1 * _tmp9;
  const Scalar _tmp15 = Scalar(0.20999999999999999) * _tmp13 - Scalar(0.20999999999999999) * _tmp14;
  const Scalar _tmp16 = _tmp12 - _tmp15;
  const Scalar _tmp17 = _tmp16 + _tmp5;
  const Scalar _tmp18 = _tmp17 + p_init0;
  const Scalar _tmp19 = Scalar(0.20999999999999999) * _tmp13 + Scalar(0.20999999999999999) * _tmp14;
  const Scalar _tmp20 = _tmp3 * _tmp7;
  const Scalar _tmp21 = _tmp6 * _tmp9;
  const Scalar _tmp22 = _tmp20 - _tmp21;
  const Scalar _tmp23 = -Scalar(0.010999999999999999) * _tmp22;
  const Scalar _tmp24 = -2 * std::pow(_tmp6, Scalar(2));
  const Scalar _tmp25 = Scalar(0.20999999999999999) * _tmp2 + Scalar(0.20999999999999999) * _tmp24 +
                        Scalar(0.20999999999999999);
  const Scalar _tmp26 = _tmp23 - _tmp25;
  const Scalar _tmp27 = _tmp19 + _tmp26;
  const Scalar _tmp28 = _tmp27 + p_init1;
  const Scalar _tmp29 = std::asinh(_tmp0 * fv1);
  const Scalar _tmp30 = Scalar(9.6622558468725703) * fh1;
  const Scalar _tmp31 = -Scalar(0.20999999999999999) * _tmp10 + Scalar(0.20999999999999999) * _tmp8;
  const Scalar _tmp32 =
      -Scalar(0.010999999999999999) * _tmp24 - Scalar(0.010999999999999999) * _tmp4;
  const Scalar _tmp33 = Scalar(0.20999999999999999) * _tmp20 + Scalar(0.20999999999999999) * _tmp21;
  const Scalar _tmp34 = _tmp32 - _tmp33;
  const Scalar _tmp35 = _tmp31 + _tmp34;
  const Scalar _tmp36 = -_tmp5;
  const Scalar _tmp37 = _tmp12 + _tmp15;
  const Scalar _tmp38 = _tmp36 + _tmp37;
  const Scalar _tmp39 = _tmp38 + p_init0;
  const Scalar _tmp40 = -_tmp19;
  const Scalar _tmp41 = _tmp23 + _tmp25;
  const Scalar _tmp42 = _tmp40 + _tmp41;
  const Scalar _tmp43 = _tmp42 + p_init1;
  const Scalar _tmp44 = Scalar(43.164000000000001) - fv1;
  const Scalar _tmp45 = Scalar(1.0) * _tmp38;
  const Scalar _tmp46 = Scalar(1.0) * _tmp42;
  const Scalar _tmp47 = -_tmp46;
  const Scalar _tmp48 = _tmp19 + _tmp41;
  const Scalar _tmp49 = Scalar(1.0) / (_tmp47 + _tmp48);
  const Scalar _tmp50 = _tmp37 + _tmp5;
  const Scalar _tmp51 = _tmp49 * (_tmp45 - _tmp50);
  const Scalar _tmp52 = _tmp45 + _tmp46 * _tmp51;
  const Scalar _tmp53 = -_tmp31;
  const Scalar _tmp54 = _tmp32 + _tmp33;
  const Scalar _tmp55 = _tmp53 + _tmp54;
  const Scalar _tmp56 = _tmp26 + _tmp40;
  const Scalar _tmp57 = _tmp56 + p_init1;
  const Scalar _tmp58 = _tmp57 + Scalar(8.3196563700000006);
  const Scalar _tmp59 = _tmp16 + _tmp36;
  const Scalar _tmp60 = _tmp59 + p_init0;
  const Scalar _tmp61 = _tmp60 + Scalar(1.9874742000000001);
  const Scalar _tmp62 = std::pow(Scalar(std::pow(_tmp58, Scalar(2)) + std::pow(_tmp61, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp63 = _tmp61 * _tmp62;
  const Scalar _tmp64 = _tmp34 + _tmp53;
  const Scalar _tmp65 = _tmp31 + _tmp54;
  const Scalar _tmp66 = _tmp50 + p_init0;
  const Scalar _tmp67 = _tmp66 + Scalar(-2.71799795);
  const Scalar _tmp68 = _tmp48 + p_init1;
  const Scalar _tmp69 = _tmp68 + Scalar(-4.7752063900000001);
  const Scalar _tmp70 = std::pow(Scalar(std::pow(_tmp67, Scalar(2)) + std::pow(_tmp69, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp71 = _tmp67 * _tmp70;
  const Scalar _tmp72 = _tmp55 * _tmp71 - _tmp65 * _tmp71;
  const Scalar _tmp73 = _tmp69 * _tmp70;
  const Scalar _tmp74 = _tmp43 + Scalar(-4.8333311099999996);
  const Scalar _tmp75 = _tmp39 + Scalar(1.79662371);
  const Scalar _tmp76 = Scalar(1.0) / (_tmp75);
  const Scalar _tmp77 = _tmp74 * _tmp76;
  const Scalar _tmp78 = Scalar(1.0) / (_tmp71 * _tmp77 - _tmp73);
  const Scalar _tmp79 = _tmp58 * _tmp62;
  const Scalar _tmp80 = _tmp63 * _tmp77 - _tmp79;
  const Scalar _tmp81 = _tmp78 * _tmp80;
  const Scalar _tmp82 = _tmp55 * _tmp77;
  const Scalar _tmp83 = _tmp78 * (_tmp65 * _tmp73 - _tmp71 * _tmp82);
  const Scalar _tmp84 = -_tmp63 * _tmp82 + _tmp64 * _tmp79 - _tmp80 * _tmp83;
  const Scalar _tmp85 = -_tmp51 * _tmp84 + _tmp55 * _tmp63 - _tmp63 * _tmp64 - _tmp72 * _tmp81;
  const Scalar _tmp86 = Scalar(1.0) / (_tmp85);
  const Scalar _tmp87 = 0;
  const Scalar _tmp88 = _tmp71 * _tmp81;
  const Scalar _tmp89 =
      std::sqrt(Scalar(std::pow(_tmp74, Scalar(2)) + std::pow(_tmp75, Scalar(2))));
  const Scalar _tmp90 = _tmp76 * _tmp89;
  const Scalar _tmp91 = _tmp77 * _tmp78;
  const Scalar _tmp92 = _tmp77 * _tmp83 + _tmp82;
  const Scalar _tmp93 = -_tmp51 * _tmp92 - _tmp55 + _tmp72 * _tmp91;
  const Scalar _tmp94 = Scalar(1.0) / (_tmp89);
  const Scalar _tmp95 = _tmp90 * (_tmp38 * _tmp74 * _tmp94 - _tmp42 * _tmp75 * _tmp94);
  const Scalar _tmp96 = _tmp48 * _tmp71 - _tmp50 * _tmp73 + _tmp71 * _tmp95;
  const Scalar _tmp97 = _tmp56 * _tmp63 - _tmp59 * _tmp79 + _tmp63 * _tmp95 - _tmp81 * _tmp96;
  const Scalar _tmp98 = _tmp86 * _tmp97;
  const Scalar _tmp99 = Scalar(1.0) / (_tmp97);
  const Scalar _tmp100 = _tmp85 * _tmp99;
  const Scalar _tmp101 = _tmp100 * (_tmp91 * _tmp96 - _tmp93 * _tmp98 - _tmp95);
  const Scalar _tmp102 = _tmp86 * (_tmp101 + _tmp93);
  const Scalar _tmp103 = -_tmp102 * _tmp80 - _tmp77;
  const Scalar _tmp104 = _tmp71 * _tmp78;
  const Scalar _tmp105 = _tmp18 + Scalar(-2.5202214700000001);
  const Scalar _tmp106 = _tmp28 + Scalar(8.3888750099999996);
  const Scalar _tmp107 =
      std::pow(Scalar(std::pow(_tmp105, Scalar(2)) + std::pow(_tmp106, Scalar(2))),
               Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp108 = _tmp105 * _tmp107;
  const Scalar _tmp109 = _tmp108 * fh1;
  const Scalar _tmp110 = Scalar(1.0) * _tmp78;
  const Scalar _tmp111 = Scalar(1.0) * _tmp83;
  const Scalar _tmp112 = -_tmp110 * _tmp72 + _tmp111 * _tmp51;
  const Scalar _tmp113 = _tmp100 * (-_tmp110 * _tmp96 - _tmp112 * _tmp98);
  const Scalar _tmp114 = _tmp86 * (_tmp112 + _tmp113);
  const Scalar _tmp115 = -_tmp114 * _tmp80 + Scalar(1.0);
  const Scalar _tmp116 = _tmp106 * _tmp107;
  const Scalar _tmp117 = _tmp116 * fh1;
  const Scalar _tmp118 = Scalar(1.0) * _tmp99;
  const Scalar _tmp119 = fh1 * (-_tmp108 * _tmp27 + _tmp116 * _tmp17);
  const Scalar _tmp120 = -_tmp109 * _tmp90 * (_tmp102 * _tmp63 + _tmp103 * _tmp104 + Scalar(1.0)) -
                         _tmp117 * _tmp90 * (_tmp104 * _tmp115 + _tmp114 * _tmp63) -
                         _tmp119 * _tmp90 * (_tmp118 * _tmp63 - _tmp118 * _tmp88) -
                         _tmp44 * _tmp90 * (_tmp63 * _tmp87 - _tmp87 * _tmp88);
  const Scalar _tmp121 = Scalar(1.0) / (_tmp120);
  const Scalar _tmp122 = _tmp47 + _tmp56;
  const Scalar _tmp123 = _tmp122 * _tmp51;
  const Scalar _tmp124 = Scalar(1.0) / (-_tmp123 + _tmp45 - _tmp59);
  const Scalar _tmp125 = Scalar(1.0) * _tmp124;
  const Scalar _tmp126 = _tmp122 * _tmp124;
  const Scalar _tmp127 = -_tmp111 + _tmp113 * _tmp126 - _tmp114 * _tmp84;
  const Scalar _tmp128 = Scalar(1.0) * _tmp49;
  const Scalar _tmp129 = _tmp101 * _tmp126 - _tmp102 * _tmp84 + _tmp92;
  const Scalar _tmp130 = _tmp125 * _tmp51;
  const Scalar _tmp131 = _tmp123 * _tmp125 + Scalar(1.0);
  const Scalar _tmp132 = _tmp35 * fh1;
  const Scalar _tmp133 = -_tmp116 * _tmp132 - Scalar(3.29616) * _tmp22 - _tmp27 * fv1;
  const Scalar _tmp134 = _tmp108 * _tmp132 + Scalar(3.29616) * _tmp11 + _tmp17 * fv1;
  const Scalar _tmp135 = _tmp122 * _tmp49;
  const Scalar _tmp136 = _tmp124 * _tmp52;
  const Scalar _tmp137 = -_tmp122 * _tmp136 + _tmp47 - _tmp84 * _tmp87;
  const Scalar _tmp138 = _tmp100 * _tmp125;
  const Scalar _tmp139 = -_tmp118 * _tmp84 + _tmp122 * _tmp138;
  const Scalar _tmp140 = std::asinh(
      _tmp121 * (Scalar(1.0) * _tmp109 * (_tmp101 * _tmp125 - _tmp128 * _tmp129) +
                 Scalar(1.0) * _tmp117 * (_tmp113 * _tmp125 - _tmp127 * _tmp128) +
                 Scalar(1.0) * _tmp119 * (-_tmp128 * _tmp139 + _tmp138) +
                 Scalar(1.0) * _tmp133 * (-_tmp128 * _tmp131 + _tmp130) +
                 Scalar(1.0) * _tmp134 * (_tmp125 * _tmp135 - _tmp125) +
                 Scalar(1.0) * _tmp44 * (-_tmp125 * _tmp52 - _tmp128 * _tmp137 + Scalar(1.0))));
  const Scalar _tmp141 = Scalar(9.6622558468725703) * _tmp120;
  const Scalar _tmp142 = _tmp118 * _tmp119;
  const Scalar _tmp143 = _tmp44 * _tmp87;
  const Scalar _tmp144 =
      _tmp103 * _tmp109 * _tmp78 + _tmp115 * _tmp117 * _tmp78 - _tmp142 * _tmp81 - _tmp143 * _tmp81;
  const Scalar _tmp145 = Scalar(1.0) / (_tmp144);
  const Scalar _tmp146 = _tmp125 * _tmp134;
  const Scalar _tmp147 =
      std::asinh(_tmp145 * (_tmp109 * _tmp129 * _tmp49 + _tmp117 * _tmp127 * _tmp49 +
                            _tmp119 * _tmp139 * _tmp49 + _tmp131 * _tmp133 * _tmp49 -
                            _tmp135 * _tmp146 + _tmp137 * _tmp44 * _tmp49));
  const Scalar _tmp148 = Scalar(9.6622558468725703) * _tmp144;
  const Scalar _tmp149 = _tmp102 * _tmp109 + _tmp114 * _tmp117 + _tmp142 + _tmp143;
  const Scalar _tmp150 = Scalar(1.0) / (_tmp149);
  const Scalar _tmp151 =
      std::asinh(_tmp150 * (-_tmp101 * _tmp109 * _tmp124 - _tmp113 * _tmp117 * _tmp124 -
                            _tmp119 * _tmp138 - _tmp130 * _tmp133 + _tmp136 * _tmp44 + _tmp146));
  const Scalar _tmp152 = Scalar(9.6622558468725703) * _tmp149;

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
      _tmp35 + p_init2;
  _res(1, 0) =
      -_tmp141 *
          (Scalar(0.86625939559540499) * _tmp121 + std::cosh(Scalar(1.0) * _tmp140) -
           std::cosh(
               Scalar(0.1034955) * _tmp121 *
               (-_tmp140 * _tmp141 -
                Scalar(4.8333311099999996) *
                    std::sqrt(Scalar(
                        std::pow(Scalar(1 - Scalar(0.20689664689659551) * _tmp43), Scalar(2)) +
                        Scalar(0.13817235445745474) *
                            std::pow(Scalar(-Scalar(0.55659957866191134) * _tmp39 - 1),
                                     Scalar(2))))))) +
      _tmp55 + p_init2;
  _res(2, 0) =
      -_tmp148 *
          (Scalar(0.86565325453551001) * _tmp145 + std::cosh(Scalar(1.0) * _tmp147) -
           std::cosh(
               Scalar(0.1034955) * _tmp145 *
               (-_tmp147 * _tmp148 -
                Scalar(4.7752063900000001) *
                    std::sqrt(Scalar(
                        Scalar(0.32397683292140877) *
                            std::pow(Scalar(1 - Scalar(0.36791786395571047) * _tmp66), Scalar(2)) +
                        std::pow(Scalar(1 - Scalar(0.20941503221602112) * _tmp68), Scalar(2))))))) +
      _tmp65 + p_init2;
  _res(3, 0) =
      -_tmp152 *
          (Scalar(0.87679799772039002) * _tmp150 + std::cosh(Scalar(1.0) * _tmp151) -
           std::cosh(
               Scalar(0.1034955) * _tmp150 *
               (-_tmp151 * _tmp152 -
                Scalar(8.3196563700000006) *
                    std::sqrt(Scalar(
                        std::pow(Scalar(-Scalar(0.12019727204189803) * _tmp57 - 1), Scalar(2)) +
                        Scalar(0.057067943376852184) *
                            std::pow(Scalar(-Scalar(0.50315118556004401) * _tmp60 - 1),
                                     Scalar(2))))))) +
      _tmp64 + p_init2;

  return _res;
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym

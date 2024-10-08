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
 * Symbolic function: IK_residual_func_cost2_Nl10
 *
 * Args:
 *     fh1: Scalar
 *     fv1: Scalar
 *     DeltaRot: Rot3
 *     position_vector: Matrix31
 *     encoder: Matrix41
 *     p_a: Matrix31
 *     p_b: Matrix31
 *     p_c: Matrix31
 *     p_d: Matrix31
 *     Rot_init: Rot3
 *     epsilon: Scalar
 *
 * Outputs:
 *     res: Matrix41
 */
template <typename Scalar>
Eigen::Matrix<Scalar, 4, 1> IkResidualFuncCost2Nl10(
    const Scalar fh1, const Scalar fv1, const sym::Rot3<Scalar>& DeltaRot,
    const Eigen::Matrix<Scalar, 3, 1>& position_vector, const Eigen::Matrix<Scalar, 4, 1>& encoder,
    const Eigen::Matrix<Scalar, 3, 1>& p_a, const Eigen::Matrix<Scalar, 3, 1>& p_b,
    const Eigen::Matrix<Scalar, 3, 1>& p_c, const Eigen::Matrix<Scalar, 3, 1>& p_d,
    const sym::Rot3<Scalar>& Rot_init, const Scalar epsilon) {
  // Total ops: 493

  // Unused inputs
  (void)encoder;
  (void)epsilon;

  // Input arrays
  const Eigen::Matrix<Scalar, 4, 1>& _DeltaRot = DeltaRot.Data();
  const Eigen::Matrix<Scalar, 4, 1>& _Rot_init = Rot_init.Data();

  // Intermediate terms (161)
  const Scalar _tmp0 = Scalar(1.0) / (fh1);
  const Scalar _tmp1 = _DeltaRot[0] * _Rot_init[3] - _DeltaRot[1] * _Rot_init[2] +
                       _DeltaRot[2] * _Rot_init[1] + _DeltaRot[3] * _Rot_init[0];
  const Scalar _tmp2 = _DeltaRot[0] * _Rot_init[2] + _DeltaRot[1] * _Rot_init[3] -
                       _DeltaRot[2] * _Rot_init[0] + _DeltaRot[3] * _Rot_init[1];
  const Scalar _tmp3 = 2 * _tmp2;
  const Scalar _tmp4 = _tmp1 * _tmp3;
  const Scalar _tmp5 = -_DeltaRot[0] * _Rot_init[1] + _DeltaRot[1] * _Rot_init[0] +
                       _DeltaRot[2] * _Rot_init[3] + _DeltaRot[3] * _Rot_init[2];
  const Scalar _tmp6 = -2 * _DeltaRot[0] * _Rot_init[0] - 2 * _DeltaRot[1] * _Rot_init[1] -
                       2 * _DeltaRot[2] * _Rot_init[2] + 2 * _DeltaRot[3] * _Rot_init[3];
  const Scalar _tmp7 = _tmp5 * _tmp6;
  const Scalar _tmp8 = Scalar(0.20999999999999999) * _tmp4 + Scalar(0.20999999999999999) * _tmp7;
  const Scalar _tmp9 = _tmp3 * _tmp5;
  const Scalar _tmp10 = _tmp1 * _tmp6;
  const Scalar _tmp11 = -_tmp10 + _tmp9;
  const Scalar _tmp12 = -Scalar(0.010999999999999999) * _tmp11;
  const Scalar _tmp13 = -2 * std::pow(_tmp1, Scalar(2));
  const Scalar _tmp14 = 1 - 2 * std::pow(_tmp5, Scalar(2));
  const Scalar _tmp15 = Scalar(0.20999999999999999) * _tmp13 + Scalar(0.20999999999999999) * _tmp14;
  const Scalar _tmp16 = _tmp12 - _tmp15;
  const Scalar _tmp17 = _tmp16 + _tmp8;
  const Scalar _tmp18 = _tmp17 + position_vector(1, 0);
  const Scalar _tmp19 = Scalar(0.20999999999999999) * _tmp4 - Scalar(0.20999999999999999) * _tmp7;
  const Scalar _tmp20 = -_tmp19;
  const Scalar _tmp21 = -2 * std::pow(_tmp2, Scalar(2));
  const Scalar _tmp22 = Scalar(0.20999999999999999) * _tmp14 + Scalar(0.20999999999999999) * _tmp21;
  const Scalar _tmp23 = 2 * _tmp1 * _tmp5;
  const Scalar _tmp24 = _tmp2 * _tmp6;
  const Scalar _tmp25 = _tmp23 + _tmp24;
  const Scalar _tmp26 = -Scalar(0.010999999999999999) * _tmp25;
  const Scalar _tmp27 = _tmp22 + _tmp26;
  const Scalar _tmp28 = _tmp20 + _tmp27;
  const Scalar _tmp29 = _tmp28 + position_vector(0, 0);
  const Scalar _tmp30 = std::pow(Scalar(-_tmp18 + p_b(1, 0)), Scalar(2)) +
                        std::pow(Scalar(-_tmp29 + p_b(0, 0)), Scalar(2));
  const Scalar _tmp31 = std::asinh(_tmp0 * fv1);
  const Scalar _tmp32 = Scalar(1.4083112389913199) * fh1;
  const Scalar _tmp33 = Scalar(0.20999999999999999) * _tmp10 + Scalar(0.20999999999999999) * _tmp9;
  const Scalar _tmp34 = -_tmp33;
  const Scalar _tmp35 = -Scalar(0.010999999999999999) * _tmp13 -
                        Scalar(0.010999999999999999) * _tmp21 + Scalar(-0.010999999999999999);
  const Scalar _tmp36 = Scalar(0.20999999999999999) * _tmp23 - Scalar(0.20999999999999999) * _tmp24;
  const Scalar _tmp37 = _tmp35 + _tmp36;
  const Scalar _tmp38 = _tmp34 + _tmp37;
  const Scalar _tmp39 = _tmp35 - _tmp36;
  const Scalar _tmp40 = _tmp33 + _tmp39;
  const Scalar _tmp41 = -_tmp8;
  const Scalar _tmp42 = _tmp12 + _tmp15;
  const Scalar _tmp43 = _tmp41 + _tmp42;
  const Scalar _tmp44 = _tmp43 + position_vector(1, 0);
  const Scalar _tmp45 = -_tmp22 + _tmp26;
  const Scalar _tmp46 = _tmp19 + _tmp45;
  const Scalar _tmp47 = _tmp46 + position_vector(0, 0);
  const Scalar _tmp48 = std::pow(Scalar(-_tmp44 + p_d(1, 0)), Scalar(2)) +
                        std::pow(Scalar(-_tmp47 + p_d(0, 0)), Scalar(2));
  const Scalar _tmp49 = Scalar(333.54000000000002) - fv1;
  const Scalar _tmp50 = _tmp33 + _tmp37;
  const Scalar _tmp51 = _tmp19 + _tmp27;
  const Scalar _tmp52 = _tmp51 + position_vector(0, 0);
  const Scalar _tmp53 = _tmp52 - p_c(0, 0);
  const Scalar _tmp54 = _tmp42 + _tmp8;
  const Scalar _tmp55 = _tmp54 + position_vector(1, 0);
  const Scalar _tmp56 = _tmp55 - p_c(1, 0);
  const Scalar _tmp57 = std::pow(Scalar(std::pow(_tmp53, Scalar(2)) + std::pow(_tmp56, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp58 = _tmp53 * _tmp57;
  const Scalar _tmp59 = _tmp34 + _tmp39;
  const Scalar _tmp60 = _tmp20 + _tmp45;
  const Scalar _tmp61 = _tmp60 + position_vector(0, 0);
  const Scalar _tmp62 = _tmp61 - p_a(0, 0);
  const Scalar _tmp63 = _tmp16 + _tmp41;
  const Scalar _tmp64 = _tmp63 + position_vector(1, 0);
  const Scalar _tmp65 = _tmp64 - p_a(1, 0);
  const Scalar _tmp66 = std::pow(Scalar(std::pow(_tmp62, Scalar(2)) + std::pow(_tmp65, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp67 = _tmp62 * _tmp66;
  const Scalar _tmp68 = _tmp40 * _tmp67;
  const Scalar _tmp69 = -_tmp59 * _tmp67 + _tmp68;
  const Scalar _tmp70 = _tmp47 - p_d(0, 0);
  const Scalar _tmp71 = Scalar(1.0) / (_tmp70);
  const Scalar _tmp72 = _tmp44 - p_d(1, 0);
  const Scalar _tmp73 = _tmp71 * _tmp72;
  const Scalar _tmp74 = _tmp65 * _tmp66;
  const Scalar _tmp75 = Scalar(1.0) / (_tmp67 * _tmp73 - _tmp74);
  const Scalar _tmp76 = _tmp56 * _tmp57;
  const Scalar _tmp77 = _tmp58 * _tmp73 - _tmp76;
  const Scalar _tmp78 = _tmp75 * _tmp77;
  const Scalar _tmp79 = _tmp40 * _tmp58;
  const Scalar _tmp80 = _tmp59 * _tmp74 - _tmp68 * _tmp73;
  const Scalar _tmp81 = _tmp50 * _tmp76 - _tmp73 * _tmp79 - _tmp78 * _tmp80;
  const Scalar _tmp82 = Scalar(1.0) * _tmp43;
  const Scalar _tmp83 = -_tmp82;
  const Scalar _tmp84 = Scalar(1.0) / (_tmp63 + _tmp83);
  const Scalar _tmp85 = Scalar(1.0) * _tmp46;
  const Scalar _tmp86 = _tmp84 * (-_tmp60 + _tmp85);
  const Scalar _tmp87 = -_tmp50 * _tmp58 - _tmp69 * _tmp78 + _tmp79 - _tmp81 * _tmp86;
  const Scalar _tmp88 = Scalar(1.0) / (_tmp87);
  const Scalar _tmp89 = _tmp82 * _tmp86 + _tmp85;
  const Scalar _tmp90 = 0;
  const Scalar _tmp91 = _tmp67 * _tmp78;
  const Scalar _tmp92 =
      std::sqrt(Scalar(std::pow(_tmp70, Scalar(2)) + std::pow(_tmp72, Scalar(2))));
  const Scalar _tmp93 = _tmp71 * _tmp92;
  const Scalar _tmp94 = Scalar(1.0) * _tmp75;
  const Scalar _tmp95 = _tmp80 * _tmp94;
  const Scalar _tmp96 = -_tmp69 * _tmp94 + _tmp86 * _tmp95;
  const Scalar _tmp97 = Scalar(1.0) / (_tmp92);
  const Scalar _tmp98 = _tmp93 * (-_tmp43 * _tmp70 * _tmp97 + _tmp46 * _tmp72 * _tmp97);
  const Scalar _tmp99 = -_tmp60 * _tmp74 + _tmp63 * _tmp67 + _tmp67 * _tmp98;
  const Scalar _tmp100 = -_tmp51 * _tmp76 + _tmp54 * _tmp58 + _tmp58 * _tmp98 - _tmp78 * _tmp99;
  const Scalar _tmp101 = _tmp100 * _tmp88;
  const Scalar _tmp102 = Scalar(1.0) / (_tmp100);
  const Scalar _tmp103 = _tmp102 * _tmp87;
  const Scalar _tmp104 = _tmp103 * (-_tmp101 * _tmp96 - _tmp94 * _tmp99);
  const Scalar _tmp105 = _tmp104 + _tmp96;
  const Scalar _tmp106 = _tmp58 * _tmp88;
  const Scalar _tmp107 = _tmp77 * _tmp88;
  const Scalar _tmp108 = -_tmp105 * _tmp107 + Scalar(1.0);
  const Scalar _tmp109 = _tmp67 * _tmp75;
  const Scalar _tmp110 = _tmp29 - p_b(0, 0);
  const Scalar _tmp111 = _tmp18 - p_b(1, 0);
  const Scalar _tmp112 =
      std::pow(Scalar(std::pow(_tmp110, Scalar(2)) + std::pow(_tmp111, Scalar(2))),
               Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp113 = _tmp111 * _tmp112;
  const Scalar _tmp114 = _tmp113 * fh1;
  const Scalar _tmp115 = _tmp73 * _tmp75;
  const Scalar _tmp116 = _tmp115 * _tmp80 + _tmp40 * _tmp73;
  const Scalar _tmp117 = _tmp115 * _tmp69 - _tmp116 * _tmp86 - _tmp40;
  const Scalar _tmp118 = _tmp103 * (-_tmp101 * _tmp117 + _tmp115 * _tmp99 - _tmp98);
  const Scalar _tmp119 = _tmp117 + _tmp118;
  const Scalar _tmp120 = -_tmp107 * _tmp119 - _tmp73;
  const Scalar _tmp121 = _tmp110 * _tmp112;
  const Scalar _tmp122 = _tmp121 * fh1;
  const Scalar _tmp123 = Scalar(1.0) * _tmp102;
  const Scalar _tmp124 = fh1 * (_tmp113 * _tmp28 - _tmp121 * _tmp17);
  const Scalar _tmp125 = -_tmp114 * _tmp93 * (_tmp105 * _tmp106 + _tmp108 * _tmp109) -
                         _tmp122 * _tmp93 * (_tmp106 * _tmp119 + _tmp109 * _tmp120 + Scalar(1.0)) -
                         _tmp124 * _tmp93 * (_tmp123 * _tmp58 - _tmp123 * _tmp91) -
                         _tmp49 * _tmp93 * (_tmp58 * _tmp90 - _tmp90 * _tmp91);
  const Scalar _tmp126 = Scalar(1.0) / (_tmp125);
  const Scalar _tmp127 = _tmp54 + _tmp83;
  const Scalar _tmp128 = _tmp127 * _tmp86;
  const Scalar _tmp129 = Scalar(1.0) / (-_tmp128 - _tmp51 + _tmp85);
  const Scalar _tmp130 = Scalar(1.0) * _tmp129;
  const Scalar _tmp131 = _tmp127 * _tmp84;
  const Scalar _tmp132 = _tmp38 * fh1;
  const Scalar _tmp133 = _tmp121 * _tmp132 + Scalar(40.024799999999999) * _tmp25 + _tmp28 * fv1;
  const Scalar _tmp134 = -Scalar(40.024799999999999) * _tmp11 - _tmp113 * _tmp132 - _tmp17 * fv1;
  const Scalar _tmp135 = _tmp130 * _tmp86;
  const Scalar _tmp136 = _tmp128 * _tmp130 + Scalar(1.0);
  const Scalar _tmp137 = Scalar(1.0) * _tmp84;
  const Scalar _tmp138 = _tmp103 * _tmp130;
  const Scalar _tmp139 = -_tmp123 * _tmp81 + _tmp127 * _tmp138;
  const Scalar _tmp140 = _tmp129 * _tmp89;
  const Scalar _tmp141 = -_tmp127 * _tmp140 - _tmp81 * _tmp90 + _tmp83;
  const Scalar _tmp142 = _tmp81 * _tmp88;
  const Scalar _tmp143 = _tmp127 * _tmp129;
  const Scalar _tmp144 = _tmp116 + _tmp118 * _tmp143 - _tmp119 * _tmp142;
  const Scalar _tmp145 = _tmp104 * _tmp143 - _tmp105 * _tmp142 - _tmp95;
  const Scalar _tmp146 = std::asinh(
      _tmp126 * (Scalar(1.0) * _tmp114 * (_tmp104 * _tmp130 - _tmp137 * _tmp145) +
                 Scalar(1.0) * _tmp122 * (_tmp118 * _tmp130 - _tmp137 * _tmp144) +
                 Scalar(1.0) * _tmp124 * (-_tmp137 * _tmp139 + _tmp138) +
                 Scalar(1.0) * _tmp133 * (_tmp130 * _tmp131 - _tmp130) +
                 Scalar(1.0) * _tmp134 * (_tmp135 - _tmp136 * _tmp137) +
                 Scalar(1.0) * _tmp49 * (-_tmp130 * _tmp89 - _tmp137 * _tmp141 + Scalar(1.0))));
  const Scalar _tmp147 = Scalar(1.4083112389913199) * _tmp125;
  const Scalar _tmp148 = std::pow(Scalar(-_tmp61 + p_a(0, 0)), Scalar(2)) +
                         std::pow(Scalar(-_tmp64 + p_a(1, 0)), Scalar(2));
  const Scalar _tmp149 = _tmp123 * _tmp124;
  const Scalar _tmp150 = _tmp49 * _tmp90;
  const Scalar _tmp151 =
      _tmp108 * _tmp114 * _tmp75 + _tmp120 * _tmp122 * _tmp75 - _tmp149 * _tmp78 - _tmp150 * _tmp78;
  const Scalar _tmp152 = Scalar(1.0) / (_tmp151);
  const Scalar _tmp153 = _tmp130 * _tmp133;
  const Scalar _tmp154 =
      std::asinh(_tmp152 * (_tmp114 * _tmp145 * _tmp84 + _tmp122 * _tmp144 * _tmp84 +
                            _tmp124 * _tmp139 * _tmp84 - _tmp131 * _tmp153 +
                            _tmp134 * _tmp136 * _tmp84 + _tmp141 * _tmp49 * _tmp84));
  const Scalar _tmp155 = Scalar(1.4083112389913199) * _tmp151;
  const Scalar _tmp156 = std::pow(Scalar(-_tmp52 + p_c(0, 0)), Scalar(2)) +
                         std::pow(Scalar(-_tmp55 + p_c(1, 0)), Scalar(2));
  const Scalar _tmp157 =
      _tmp105 * _tmp114 * _tmp88 + _tmp119 * _tmp122 * _tmp88 + _tmp149 + _tmp150;
  const Scalar _tmp158 = Scalar(1.0) / (_tmp157);
  const Scalar _tmp159 =
      std::asinh(_tmp158 * (-_tmp104 * _tmp114 * _tmp129 - _tmp118 * _tmp122 * _tmp129 -
                            _tmp124 * _tmp138 - _tmp134 * _tmp135 + _tmp140 * _tmp49 + _tmp153));
  const Scalar _tmp160 = Scalar(1.4083112389913199) * _tmp157;

  // Output terms (1)
  Eigen::Matrix<Scalar, 4, 1> _res;

  _res(0, 0) =
      _tmp32 *
          (-std::sinh(Scalar(1.0) * _tmp31) - std::sinh(Scalar(0.71007031138673404) * _tmp0 *
                                                        (-std::sqrt(_tmp30) - _tmp31 * _tmp32))) -
      std::sqrt(Scalar(_tmp30 +
                       std::pow(Scalar(-_tmp38 + p_b(2, 0) - position_vector(2, 0)), Scalar(2))));
  _res(1, 0) =
      _tmp147 * (-std::sinh(Scalar(1.0) * _tmp146) -
                 std::sinh(Scalar(0.71007031138673404) * _tmp126 *
                           (-_tmp146 * _tmp147 - std::sqrt(_tmp48)))) -
      std::sqrt(Scalar(_tmp48 +
                       std::pow(Scalar(-_tmp40 + p_d(2, 0) - position_vector(2, 0)), Scalar(2))));
  _res(2, 0) =
      _tmp155 * (-std::sinh(Scalar(1.0) * _tmp154) -
                 std::sinh(Scalar(0.71007031138673404) * _tmp152 *
                           (-std::sqrt(_tmp148) - _tmp154 * _tmp155))) -
      std::sqrt(Scalar(_tmp148 +
                       std::pow(Scalar(-_tmp59 + p_a(2, 0) - position_vector(2, 0)), Scalar(2))));
  _res(3, 0) =
      _tmp160 * (-std::sinh(Scalar(1.0) * _tmp159) -
                 std::sinh(Scalar(0.71007031138673404) * _tmp158 *
                           (-std::sqrt(_tmp156) - _tmp159 * _tmp160))) -
      std::sqrt(Scalar(_tmp156 +
                       std::pow(Scalar(-_tmp50 + p_c(2, 0) - position_vector(2, 0)), Scalar(2))));

  return _res;
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym

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
 * Symbolic function: IK_residual_func_cost2_Nl6
 *
 * Args:
 *     fh1: Scalar
 *     fv1: Scalar
 *     DeltaRot: Rot3
 *     position_vector: Matrix31
 *     Rot_init: Rot3
 *     epsilon: Scalar
 *
 * Outputs:
 *     res: Matrix41
 */
template <typename Scalar>
Eigen::Matrix<Scalar, 4, 1> IkResidualFuncCost2Nl6(
    const Scalar fh1, const Scalar fv1, const sym::Rot3<Scalar>& DeltaRot,
    const Eigen::Matrix<Scalar, 3, 1>& position_vector, const sym::Rot3<Scalar>& Rot_init,
    const Scalar epsilon) {
  // Total ops: 530

  // Unused inputs
  (void)epsilon;

  // Input arrays
  const Eigen::Matrix<Scalar, 4, 1>& _DeltaRot = DeltaRot.Data();
  const Eigen::Matrix<Scalar, 4, 1>& _Rot_init = Rot_init.Data();

  // Intermediate terms (162)
  const Scalar _tmp0 = Scalar(1.0) / (fh1);
  const Scalar _tmp1 = _DeltaRot[0] * _Rot_init[3] - _DeltaRot[1] * _Rot_init[2] +
                       _DeltaRot[2] * _Rot_init[1] + _DeltaRot[3] * _Rot_init[0];
  const Scalar _tmp2 = -2 * std::pow(_tmp1, Scalar(2));
  const Scalar _tmp3 = -_DeltaRot[0] * _Rot_init[1] + _DeltaRot[1] * _Rot_init[0] +
                       _DeltaRot[2] * _Rot_init[3] + _DeltaRot[3] * _Rot_init[2];
  const Scalar _tmp4 = -2 * std::pow(_tmp3, Scalar(2));
  const Scalar _tmp5 = Scalar(0.20999999999999999) * _tmp2 + Scalar(0.20999999999999999) * _tmp4 +
                       Scalar(0.20999999999999999);
  const Scalar _tmp6 = -_tmp5;
  const Scalar _tmp7 = _DeltaRot[0] * _Rot_init[2] + _DeltaRot[1] * _Rot_init[3] -
                       _DeltaRot[2] * _Rot_init[0] + _DeltaRot[3] * _Rot_init[1];
  const Scalar _tmp8 = 2 * _tmp3 * _tmp7;
  const Scalar _tmp9 = -2 * _DeltaRot[0] * _Rot_init[0] - 2 * _DeltaRot[1] * _Rot_init[1] -
                       2 * _DeltaRot[2] * _Rot_init[2] + 2 * _DeltaRot[3] * _Rot_init[3];
  const Scalar _tmp10 = _tmp1 * _tmp9;
  const Scalar _tmp11 = -_tmp10 + _tmp8;
  const Scalar _tmp12 = -Scalar(0.010999999999999999) * _tmp11;
  const Scalar _tmp13 = 2 * _tmp1;
  const Scalar _tmp14 = _tmp13 * _tmp7;
  const Scalar _tmp15 = _tmp3 * _tmp9;
  const Scalar _tmp16 = Scalar(0.20999999999999999) * _tmp14 + Scalar(0.20999999999999999) * _tmp15;
  const Scalar _tmp17 = _tmp12 + _tmp16;
  const Scalar _tmp18 = _tmp17 + _tmp6;
  const Scalar _tmp19 = _tmp18 + position_vector(1, 0);
  const Scalar _tmp20 = 1 - 2 * std::pow(_tmp7, Scalar(2));
  const Scalar _tmp21 = Scalar(0.20999999999999999) * _tmp20 + Scalar(0.20999999999999999) * _tmp4;
  const Scalar _tmp22 = _tmp13 * _tmp3;
  const Scalar _tmp23 = _tmp7 * _tmp9;
  const Scalar _tmp24 = _tmp22 + _tmp23;
  const Scalar _tmp25 = -Scalar(0.010999999999999999) * _tmp24;
  const Scalar _tmp26 = Scalar(0.20999999999999999) * _tmp14 - Scalar(0.20999999999999999) * _tmp15;
  const Scalar _tmp27 = _tmp25 - _tmp26;
  const Scalar _tmp28 = _tmp21 + _tmp27;
  const Scalar _tmp29 = _tmp28 + position_vector(0, 0);
  const Scalar _tmp30 =
      Scalar(15625.0) * std::pow(Scalar(1 - Scalar(0.0080000000000000002) * _tmp29), Scalar(2)) +
      Scalar(12100.0) * std::pow(Scalar(-Scalar(0.0090909090909090905) * _tmp19 - 1), Scalar(2));
  const Scalar _tmp31 = std::asinh(_tmp0 * fv1);
  const Scalar _tmp32 = Scalar(1.4083112389913199) * fh1;
  const Scalar _tmp33 = Scalar(0.20999999999999999) * _tmp22 - Scalar(0.20999999999999999) * _tmp23;
  const Scalar _tmp34 =
      -Scalar(0.010999999999999999) * _tmp2 - Scalar(0.010999999999999999) * _tmp20;
  const Scalar _tmp35 = Scalar(0.20999999999999999) * _tmp10 + Scalar(0.20999999999999999) * _tmp8;
  const Scalar _tmp36 = _tmp34 - _tmp35;
  const Scalar _tmp37 = _tmp33 + _tmp36;
  const Scalar _tmp38 = -_tmp21;
  const Scalar _tmp39 = _tmp27 + _tmp38;
  const Scalar _tmp40 = Scalar(1.0) * _tmp39;
  const Scalar _tmp41 = _tmp25 + _tmp26;
  const Scalar _tmp42 = _tmp21 + _tmp41;
  const Scalar _tmp43 = _tmp12 - _tmp16;
  const Scalar _tmp44 = _tmp43 + _tmp6;
  const Scalar _tmp45 = Scalar(1.0) * _tmp44;
  const Scalar _tmp46 = -_tmp45;
  const Scalar _tmp47 = _tmp17 + _tmp5;
  const Scalar _tmp48 = Scalar(1.0) / (_tmp46 + _tmp47);
  const Scalar _tmp49 = _tmp48 * (_tmp40 - _tmp42);
  const Scalar _tmp50 = _tmp34 + _tmp35;
  const Scalar _tmp51 = _tmp33 + _tmp50;
  const Scalar _tmp52 = _tmp47 + position_vector(1, 0);
  const Scalar _tmp53 = _tmp52 + Scalar(-110.0);
  const Scalar _tmp54 = _tmp42 + position_vector(0, 0);
  const Scalar _tmp55 = _tmp54 + Scalar(-125.0);
  const Scalar _tmp56 = std::pow(Scalar(std::pow(_tmp53, Scalar(2)) + std::pow(_tmp55, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp57 = _tmp53 * _tmp56;
  const Scalar _tmp58 = -_tmp33;
  const Scalar _tmp59 = _tmp36 + _tmp58;
  const Scalar _tmp60 = _tmp44 + position_vector(1, 0);
  const Scalar _tmp61 = _tmp60 + Scalar(110.0);
  const Scalar _tmp62 = _tmp39 + position_vector(0, 0);
  const Scalar _tmp63 = _tmp62 + Scalar(125.0);
  const Scalar _tmp64 = Scalar(1.0) / (_tmp63);
  const Scalar _tmp65 = _tmp61 * _tmp64;
  const Scalar _tmp66 = _tmp59 * _tmp65;
  const Scalar _tmp67 = _tmp55 * _tmp56;
  const Scalar _tmp68 = _tmp51 * _tmp57 - _tmp66 * _tmp67;
  const Scalar _tmp69 = Scalar(1.0) / (-_tmp57 + _tmp65 * _tmp67);
  const Scalar _tmp70 = Scalar(1.0) * _tmp69;
  const Scalar _tmp71 = _tmp68 * _tmp70;
  const Scalar _tmp72 = -_tmp51 * _tmp67 + _tmp59 * _tmp67;
  const Scalar _tmp73 = _tmp49 * _tmp71 - _tmp70 * _tmp72;
  const Scalar _tmp74 = _tmp50 + _tmp58;
  const Scalar _tmp75 = _tmp38 + _tmp41;
  const Scalar _tmp76 = _tmp75 + position_vector(0, 0);
  const Scalar _tmp77 = _tmp76 + Scalar(125.0);
  const Scalar _tmp78 = _tmp43 + _tmp5;
  const Scalar _tmp79 = _tmp78 + position_vector(1, 0);
  const Scalar _tmp80 = _tmp79 + Scalar(-110.0);
  const Scalar _tmp81 = std::pow(Scalar(std::pow(_tmp77, Scalar(2)) + std::pow(_tmp80, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp82 = _tmp77 * _tmp81;
  const Scalar _tmp83 = _tmp59 * _tmp82;
  const Scalar _tmp84 = _tmp80 * _tmp81;
  const Scalar _tmp85 = _tmp65 * _tmp82 - _tmp84;
  const Scalar _tmp86 = _tmp69 * _tmp85;
  const Scalar _tmp87 = -_tmp65 * _tmp83 - _tmp68 * _tmp86 + _tmp74 * _tmp84;
  const Scalar _tmp88 = -_tmp49 * _tmp87 - _tmp72 * _tmp86 - _tmp74 * _tmp82 + _tmp83;
  const Scalar _tmp89 = Scalar(1.0) / (_tmp88);
  const Scalar _tmp90 =
      std::sqrt(Scalar(std::pow(_tmp61, Scalar(2)) + std::pow(_tmp63, Scalar(2))));
  const Scalar _tmp91 = Scalar(1.0) / (_tmp90);
  const Scalar _tmp92 = _tmp64 * _tmp90;
  const Scalar _tmp93 = _tmp92 * (_tmp39 * _tmp61 * _tmp91 - _tmp44 * _tmp63 * _tmp91);
  const Scalar _tmp94 = -_tmp42 * _tmp57 + _tmp47 * _tmp67 + _tmp67 * _tmp93;
  const Scalar _tmp95 = -_tmp75 * _tmp84 + _tmp78 * _tmp82 + _tmp82 * _tmp93 - _tmp86 * _tmp94;
  const Scalar _tmp96 = _tmp89 * _tmp95;
  const Scalar _tmp97 = Scalar(1.0) / (_tmp95);
  const Scalar _tmp98 = _tmp88 * _tmp97;
  const Scalar _tmp99 = _tmp98 * (-_tmp70 * _tmp94 - _tmp73 * _tmp96);
  const Scalar _tmp100 = _tmp73 + _tmp99;
  const Scalar _tmp101 = _tmp82 * _tmp89;
  const Scalar _tmp102 = _tmp85 * _tmp89;
  const Scalar _tmp103 = -_tmp100 * _tmp102 + Scalar(1.0);
  const Scalar _tmp104 = _tmp67 * _tmp69;
  const Scalar _tmp105 = _tmp19 + Scalar(110.0);
  const Scalar _tmp106 = _tmp29 + Scalar(-125.0);
  const Scalar _tmp107 =
      std::pow(Scalar(std::pow(_tmp105, Scalar(2)) + std::pow(_tmp106, Scalar(2))),
               Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp108 = _tmp105 * _tmp107;
  const Scalar _tmp109 = _tmp108 * fh1;
  const Scalar _tmp110 = Scalar(333.54000000000002) - fv1;
  const Scalar _tmp111 = _tmp40 + _tmp45 * _tmp49;
  const Scalar _tmp112 = 0;
  const Scalar _tmp113 = _tmp112 * _tmp89;
  const Scalar _tmp114 = _tmp67 * _tmp86;
  const Scalar _tmp115 = Scalar(1.0) * _tmp97;
  const Scalar _tmp116 = _tmp106 * _tmp107;
  const Scalar _tmp117 = fh1 * (_tmp108 * _tmp28 - _tmp116 * _tmp18);
  const Scalar _tmp118 = _tmp65 * _tmp69;
  const Scalar _tmp119 = _tmp118 * _tmp68 + _tmp66;
  const Scalar _tmp120 = _tmp118 * _tmp72 - _tmp119 * _tmp49 - _tmp59;
  const Scalar _tmp121 = _tmp98 * (_tmp118 * _tmp94 - _tmp120 * _tmp96 - _tmp93);
  const Scalar _tmp122 = _tmp120 + _tmp121;
  const Scalar _tmp123 = -_tmp102 * _tmp122 - _tmp65;
  const Scalar _tmp124 = _tmp116 * fh1;
  const Scalar _tmp125 = -_tmp109 * _tmp92 * (_tmp100 * _tmp101 + _tmp103 * _tmp104) -
                         _tmp110 * _tmp92 * (-_tmp113 * _tmp114 + _tmp113 * _tmp82) -
                         _tmp117 * _tmp92 * (-_tmp114 * _tmp115 + _tmp115 * _tmp82) -
                         _tmp124 * _tmp92 * (_tmp101 * _tmp122 + _tmp104 * _tmp123 + Scalar(1.0));
  const Scalar _tmp126 = Scalar(1.0) / (_tmp125);
  const Scalar _tmp127 =
      Scalar(12100.0) * std::pow(Scalar(-Scalar(0.0090909090909090905) * _tmp60 - 1), Scalar(2)) +
      Scalar(15625.0) * std::pow(Scalar(-Scalar(0.0080000000000000002) * _tmp62 - 1), Scalar(2));
  const Scalar _tmp128 = _tmp46 + _tmp78;
  const Scalar _tmp129 = _tmp128 * _tmp49;
  const Scalar _tmp130 = Scalar(1.0) / (-_tmp129 + _tmp40 - _tmp75);
  const Scalar _tmp131 = Scalar(1.0) * _tmp130;
  const Scalar _tmp132 = _tmp48 * (_tmp129 * _tmp131 + Scalar(1.0));
  const Scalar _tmp133 = _tmp131 * _tmp49;
  const Scalar _tmp134 = _tmp37 * fh1;
  const Scalar _tmp135 = -_tmp108 * _tmp134 - Scalar(40.024799999999999) * _tmp11 - _tmp18 * fv1;
  const Scalar _tmp136 = _tmp128 * _tmp130;
  const Scalar _tmp137 = _tmp87 * _tmp89;
  const Scalar _tmp138 = _tmp119 + _tmp121 * _tmp136 - _tmp122 * _tmp137;
  const Scalar _tmp139 = Scalar(1.0) * _tmp48;
  const Scalar _tmp140 = _tmp128 * _tmp48;
  const Scalar _tmp141 = _tmp116 * _tmp134 + Scalar(40.024799999999999) * _tmp24 + _tmp28 * fv1;
  const Scalar _tmp142 = _tmp111 * _tmp130;
  const Scalar _tmp143 = -_tmp112 * _tmp137 - _tmp128 * _tmp142 + _tmp46;
  const Scalar _tmp144 = _tmp131 * _tmp98;
  const Scalar _tmp145 = -_tmp115 * _tmp87 + _tmp128 * _tmp144;
  const Scalar _tmp146 = -_tmp100 * _tmp137 + _tmp136 * _tmp99 - _tmp71;
  const Scalar _tmp147 = std::asinh(
      _tmp126 * (Scalar(1.0) * _tmp109 * (_tmp131 * _tmp99 - _tmp139 * _tmp146) +
                 Scalar(1.0) * _tmp110 * (-_tmp111 * _tmp131 - _tmp139 * _tmp143 + Scalar(1.0)) +
                 Scalar(1.0) * _tmp117 * (-_tmp139 * _tmp145 + _tmp144) +
                 Scalar(1.0) * _tmp124 * (_tmp121 * _tmp131 - _tmp138 * _tmp139) +
                 Scalar(1.0) * _tmp135 * (-Scalar(1.0) * _tmp132 + _tmp133) +
                 Scalar(1.0) * _tmp141 * (_tmp131 * _tmp140 - _tmp131)));
  const Scalar _tmp148 = Scalar(1.4083112389913199) * _tmp125;
  const Scalar _tmp149 = _tmp131 * _tmp141;
  const Scalar _tmp150 = _tmp115 * _tmp117;
  const Scalar _tmp151 = _tmp110 * _tmp113;
  const Scalar _tmp152 =
      _tmp103 * _tmp109 * _tmp69 + _tmp123 * _tmp124 * _tmp69 - _tmp150 * _tmp86 - _tmp151 * _tmp86;
  const Scalar _tmp153 = Scalar(1.0) / (_tmp152);
  const Scalar _tmp154 =
      std::asinh(_tmp153 * (_tmp109 * _tmp146 * _tmp48 + _tmp110 * _tmp143 * _tmp48 +
                            _tmp117 * _tmp145 * _tmp48 + _tmp124 * _tmp138 * _tmp48 +
                            _tmp132 * _tmp135 - _tmp140 * _tmp149));
  const Scalar _tmp155 =
      Scalar(12100.0) * std::pow(Scalar(1 - Scalar(0.0090909090909090905) * _tmp52), Scalar(2)) +
      Scalar(15625.0) * std::pow(Scalar(1 - Scalar(0.0080000000000000002) * _tmp54), Scalar(2));
  const Scalar _tmp156 = Scalar(1.4083112389913199) * _tmp152;
  const Scalar _tmp157 =
      Scalar(12100.0) * std::pow(Scalar(1 - Scalar(0.0090909090909090905) * _tmp79), Scalar(2)) +
      Scalar(15625.0) * std::pow(Scalar(-Scalar(0.0080000000000000002) * _tmp76 - 1), Scalar(2));
  const Scalar _tmp158 =
      _tmp100 * _tmp109 * _tmp89 + _tmp122 * _tmp124 * _tmp89 + _tmp150 + _tmp151;
  const Scalar _tmp159 = Scalar(1.0) / (_tmp158);
  const Scalar _tmp160 =
      std::asinh(_tmp159 * (-_tmp109 * _tmp130 * _tmp99 + _tmp110 * _tmp142 - _tmp117 * _tmp144 -
                            _tmp121 * _tmp124 * _tmp130 - _tmp133 * _tmp135 + _tmp149));
  const Scalar _tmp161 = Scalar(1.4083112389913199) * _tmp158;

  // Output terms (1)
  Eigen::Matrix<Scalar, 4, 1> _res;

  _res(0, 0) =
      _tmp32 *
          (-std::sinh(Scalar(1.0) * _tmp31) - std::sinh(Scalar(0.71007031138673404) * _tmp0 *
                                                        (-std::sqrt(_tmp30) - _tmp31 * _tmp32))) -
      Scalar(48.0) * std::sqrt(Scalar(
                         Scalar(0.00043402777777777775) * _tmp30 +
                         std::pow(Scalar(-Scalar(0.020833333333333332) * _tmp37 -
                                         Scalar(0.020833333333333332) * position_vector(2, 0) + 1),
                                  Scalar(2))));
  _res(1, 0) =
      _tmp148 * (-std::sinh(Scalar(1.0) * _tmp147) -
                 std::sinh(Scalar(0.71007031138673404) * _tmp126 *
                           (-std::sqrt(_tmp127) - _tmp147 * _tmp148))) -
      Scalar(48.0) * std::sqrt(Scalar(
                         Scalar(0.00043402777777777775) * _tmp127 +
                         std::pow(Scalar(-Scalar(0.020833333333333332) * _tmp59 -
                                         Scalar(0.020833333333333332) * position_vector(2, 0) + 1),
                                  Scalar(2))));
  _res(2, 0) =
      _tmp156 * (-std::sinh(Scalar(1.0) * _tmp154) -
                 std::sinh(Scalar(0.71007031138673404) * _tmp153 *
                           (-_tmp154 * _tmp156 - std::sqrt(_tmp155)))) -
      Scalar(48.0) * std::sqrt(Scalar(
                         Scalar(0.00043402777777777775) * _tmp155 +
                         std::pow(Scalar(-Scalar(0.020833333333333332) * _tmp51 -
                                         Scalar(0.020833333333333332) * position_vector(2, 0) + 1),
                                  Scalar(2))));
  _res(3, 0) =
      _tmp161 * (-std::sinh(Scalar(1.0) * _tmp160) -
                 std::sinh(Scalar(0.71007031138673404) * _tmp159 *
                           (-std::sqrt(_tmp157) - _tmp160 * _tmp161))) -
      Scalar(48.0) * std::sqrt(Scalar(
                         Scalar(0.00043402777777777775) * _tmp157 +
                         std::pow(Scalar(-Scalar(0.020833333333333332) * _tmp74 -
                                         Scalar(0.020833333333333332) * position_vector(2, 0) + 1),
                                  Scalar(2))));

  return _res;
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym

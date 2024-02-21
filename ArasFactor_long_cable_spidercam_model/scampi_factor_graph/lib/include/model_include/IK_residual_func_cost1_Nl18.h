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
 * Symbolic function: IK_residual_func_cost1_Nl18
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
Eigen::Matrix<Scalar, 4, 1> IkResidualFuncCost1Nl18(
    const Scalar fh1, const Scalar fv1, const sym::Rot3<Scalar>& DeltaRot,
    const Eigen::Matrix<Scalar, 3, 1>& position_vector, const sym::Rot3<Scalar>& Rot_init,
    const Scalar epsilon) {
  // Total ops: 500

  // Unused inputs
  (void)epsilon;

  // Input arrays
  const Eigen::Matrix<Scalar, 4, 1>& _DeltaRot = DeltaRot.Data();
  const Eigen::Matrix<Scalar, 4, 1>& _Rot_init = Rot_init.Data();

  // Intermediate terms (153)
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
  const Scalar _tmp8 = Scalar(0.20999999999999999) * _tmp4 - Scalar(0.20999999999999999) * _tmp7;
  const Scalar _tmp9 = -2 * std::pow(_tmp5, Scalar(2));
  const Scalar _tmp10 = 1 - 2 * std::pow(_tmp2, Scalar(2));
  const Scalar _tmp11 = Scalar(0.20999999999999999) * _tmp10 + Scalar(0.20999999999999999) * _tmp9;
  const Scalar _tmp12 = 2 * _tmp1 * _tmp5;
  const Scalar _tmp13 = _tmp2 * _tmp6;
  const Scalar _tmp14 = _tmp12 + _tmp13;
  const Scalar _tmp15 = -Scalar(0.010999999999999999) * _tmp14;
  const Scalar _tmp16 = -_tmp11 + _tmp15;
  const Scalar _tmp17 = _tmp16 + _tmp8;
  const Scalar _tmp18 = _tmp17 + position_vector(0, 0);
  const Scalar _tmp19 = -2 * std::pow(_tmp1, Scalar(2));
  const Scalar _tmp20 = Scalar(0.20999999999999999) * _tmp19 + Scalar(0.20999999999999999) * _tmp9 +
                        Scalar(0.20999999999999999);
  const Scalar _tmp21 = Scalar(0.20999999999999999) * _tmp4 + Scalar(0.20999999999999999) * _tmp7;
  const Scalar _tmp22 = _tmp3 * _tmp5;
  const Scalar _tmp23 = _tmp1 * _tmp6;
  const Scalar _tmp24 = _tmp22 - _tmp23;
  const Scalar _tmp25 = -Scalar(0.010999999999999999) * _tmp24;
  const Scalar _tmp26 = -_tmp21 + _tmp25;
  const Scalar _tmp27 = _tmp20 + _tmp26;
  const Scalar _tmp28 = _tmp27 + position_vector(1, 0);
  const Scalar _tmp29 = std::asinh(_tmp0 * fv1);
  const Scalar _tmp30 = Scalar(1.4083112389913199) * fh1;
  const Scalar _tmp31 = Scalar(0.20999999999999999) * _tmp22 + Scalar(0.20999999999999999) * _tmp23;
  const Scalar _tmp32 =
      -Scalar(0.010999999999999999) * _tmp10 - Scalar(0.010999999999999999) * _tmp19;
  const Scalar _tmp33 = Scalar(0.20999999999999999) * _tmp12 - Scalar(0.20999999999999999) * _tmp13;
  const Scalar _tmp34 = _tmp32 - _tmp33;
  const Scalar _tmp35 = _tmp31 + _tmp34;
  const Scalar _tmp36 = -_tmp20;
  const Scalar _tmp37 = _tmp26 + _tmp36;
  const Scalar _tmp38 = _tmp37 + position_vector(1, 0);
  const Scalar _tmp39 = _tmp38 + Scalar(110.0);
  const Scalar _tmp40 = -_tmp8;
  const Scalar _tmp41 = _tmp16 + _tmp40;
  const Scalar _tmp42 = _tmp41 + position_vector(0, 0);
  const Scalar _tmp43 = _tmp42 + Scalar(125.0);
  const Scalar _tmp44 = Scalar(1.0) / (_tmp43);
  const Scalar _tmp45 = _tmp39 * _tmp44;
  const Scalar _tmp46 = _tmp21 + _tmp25;
  const Scalar _tmp47 = _tmp20 + _tmp46;
  const Scalar _tmp48 = _tmp47 + position_vector(1, 0);
  const Scalar _tmp49 = _tmp48 + Scalar(-110.0);
  const Scalar _tmp50 = _tmp11 + _tmp15;
  const Scalar _tmp51 = _tmp50 + _tmp8;
  const Scalar _tmp52 = _tmp51 + position_vector(0, 0);
  const Scalar _tmp53 = _tmp52 + Scalar(-125.0);
  const Scalar _tmp54 = std::pow(Scalar(std::pow(_tmp49, Scalar(2)) + std::pow(_tmp53, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp55 = _tmp53 * _tmp54;
  const Scalar _tmp56 = _tmp49 * _tmp54;
  const Scalar _tmp57 = _tmp45 * _tmp55 - _tmp56;
  const Scalar _tmp58 = -_tmp31;
  const Scalar _tmp59 = _tmp34 + _tmp58;
  const Scalar _tmp60 = _tmp55 * _tmp59;
  const Scalar _tmp61 = _tmp32 + _tmp33;
  const Scalar _tmp62 = _tmp31 + _tmp61;
  const Scalar _tmp63 = _tmp58 + _tmp61;
  const Scalar _tmp64 = _tmp36 + _tmp46;
  const Scalar _tmp65 = _tmp64 + position_vector(1, 0);
  const Scalar _tmp66 = _tmp65 + Scalar(110.0);
  const Scalar _tmp67 = _tmp40 + _tmp50;
  const Scalar _tmp68 = _tmp67 + position_vector(0, 0);
  const Scalar _tmp69 = _tmp68 + Scalar(-125.0);
  const Scalar _tmp70 = std::pow(Scalar(std::pow(_tmp66, Scalar(2)) + std::pow(_tmp69, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp71 = _tmp66 * _tmp70;
  const Scalar _tmp72 = _tmp69 * _tmp70;
  const Scalar _tmp73 = _tmp59 * _tmp72;
  const Scalar _tmp74 = -_tmp45 * _tmp73 + _tmp63 * _tmp71;
  const Scalar _tmp75 = Scalar(1.0) / (_tmp45 * _tmp72 - _tmp71);
  const Scalar _tmp76 = _tmp57 * _tmp75;
  const Scalar _tmp77 = -_tmp45 * _tmp60 + _tmp56 * _tmp62 - _tmp74 * _tmp76;
  const Scalar _tmp78 = Scalar(1.0) * _tmp37;
  const Scalar _tmp79 = -_tmp78;
  const Scalar _tmp80 = Scalar(1.0) / (_tmp64 + _tmp79);
  const Scalar _tmp81 = Scalar(1.0) * _tmp41;
  const Scalar _tmp82 = _tmp80 * (-_tmp67 + _tmp81);
  const Scalar _tmp83 = -_tmp63 * _tmp72 + _tmp73;
  const Scalar _tmp84 = -_tmp55 * _tmp62 + _tmp60 - _tmp76 * _tmp83 - _tmp77 * _tmp82;
  const Scalar _tmp85 = Scalar(1.0) / (_tmp84);
  const Scalar _tmp86 =
      std::sqrt(Scalar(std::pow(_tmp39, Scalar(2)) + std::pow(_tmp43, Scalar(2))));
  const Scalar _tmp87 = Scalar(1.0) / (_tmp86);
  const Scalar _tmp88 = _tmp44 * _tmp86;
  const Scalar _tmp89 = _tmp88 * (-_tmp37 * _tmp43 * _tmp87 + _tmp39 * _tmp41 * _tmp87);
  const Scalar _tmp90 = _tmp64 * _tmp72 - _tmp67 * _tmp71 + _tmp72 * _tmp89;
  const Scalar _tmp91 = _tmp45 * _tmp75;
  const Scalar _tmp92 = _tmp45 * _tmp59 + _tmp74 * _tmp91;
  const Scalar _tmp93 = -_tmp59 - _tmp82 * _tmp92 + _tmp83 * _tmp91;
  const Scalar _tmp94 = _tmp47 * _tmp55 - _tmp51 * _tmp56 + _tmp55 * _tmp89 - _tmp76 * _tmp90;
  const Scalar _tmp95 = _tmp85 * _tmp94;
  const Scalar _tmp96 = Scalar(1.0) / (_tmp94);
  const Scalar _tmp97 = _tmp84 * _tmp96;
  const Scalar _tmp98 = _tmp97 * (-_tmp89 + _tmp90 * _tmp91 - _tmp93 * _tmp95);
  const Scalar _tmp99 = _tmp85 * (_tmp93 + _tmp98);
  const Scalar _tmp100 = -_tmp45 - _tmp57 * _tmp99;
  const Scalar _tmp101 = _tmp72 * _tmp75;
  const Scalar _tmp102 = _tmp18 + Scalar(125.0);
  const Scalar _tmp103 = _tmp28 + Scalar(-110.0);
  const Scalar _tmp104 =
      std::pow(Scalar(std::pow(_tmp102, Scalar(2)) + std::pow(_tmp103, Scalar(2))),
               Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp105 = _tmp102 * _tmp104;
  const Scalar _tmp106 = _tmp105 * fh1;
  const Scalar _tmp107 = Scalar(333.54000000000002) - fv1;
  const Scalar _tmp108 = _tmp78 * _tmp82 + _tmp81;
  const Scalar _tmp109 = 0;
  const Scalar _tmp110 = Scalar(1.0) * _tmp96;
  const Scalar _tmp111 = Scalar(1.0) * _tmp75;
  const Scalar _tmp112 = _tmp111 * _tmp57 * _tmp96;
  const Scalar _tmp113 = _tmp103 * _tmp104;
  const Scalar _tmp114 = fh1 * (-_tmp105 * _tmp27 + _tmp113 * _tmp17);
  const Scalar _tmp115 = _tmp111 * _tmp74;
  const Scalar _tmp116 = -_tmp111 * _tmp83 + _tmp115 * _tmp82;
  const Scalar _tmp117 = _tmp97 * (-_tmp111 * _tmp90 - _tmp116 * _tmp95);
  const Scalar _tmp118 = _tmp85 * (_tmp116 + _tmp117);
  const Scalar _tmp119 = -_tmp118 * _tmp57 + Scalar(1.0);
  const Scalar _tmp120 = _tmp113 * fh1;
  const Scalar _tmp121 = -_tmp106 * _tmp88 * (_tmp100 * _tmp101 + _tmp55 * _tmp99 + Scalar(1.0)) -
                         _tmp107 * _tmp88 * (_tmp109 * _tmp55 - _tmp109 * _tmp72 * _tmp76) -
                         _tmp114 * _tmp88 * (_tmp110 * _tmp55 - _tmp112 * _tmp72) -
                         _tmp120 * _tmp88 * (_tmp101 * _tmp119 + _tmp118 * _tmp55);
  const Scalar _tmp122 = Scalar(1.0) / (_tmp121);
  const Scalar _tmp123 = _tmp47 + _tmp79;
  const Scalar _tmp124 = _tmp123 * _tmp82;
  const Scalar _tmp125 = Scalar(1.0) / (-_tmp124 - _tmp51 + _tmp81);
  const Scalar _tmp126 = _tmp108 * _tmp125;
  const Scalar _tmp127 = -_tmp109 * _tmp77 - _tmp123 * _tmp126 + _tmp79;
  const Scalar _tmp128 = Scalar(1.0) * _tmp80;
  const Scalar _tmp129 = Scalar(1.0) * _tmp125;
  const Scalar _tmp130 = _tmp35 * fh1;
  const Scalar _tmp131 = -_tmp113 * _tmp130 - Scalar(40.024799999999999) * _tmp24 - _tmp27 * fv1;
  const Scalar _tmp132 = _tmp124 * _tmp129 + Scalar(1.0);
  const Scalar _tmp133 = _tmp129 * _tmp82;
  const Scalar _tmp134 = _tmp105 * _tmp130 + Scalar(40.024799999999999) * _tmp14 + _tmp17 * fv1;
  const Scalar _tmp135 = _tmp123 * _tmp80;
  const Scalar _tmp136 = _tmp123 * _tmp125;
  const Scalar _tmp137 = _tmp136 * _tmp98 - _tmp77 * _tmp99 + _tmp92;
  const Scalar _tmp138 = -_tmp115 + _tmp117 * _tmp136 - _tmp118 * _tmp77;
  const Scalar _tmp139 = _tmp129 * _tmp97;
  const Scalar _tmp140 = -_tmp110 * _tmp77 + _tmp123 * _tmp139;
  const Scalar _tmp141 = std::asinh(
      _tmp122 * (Scalar(1.0) * _tmp106 * (-_tmp128 * _tmp137 + _tmp129 * _tmp98) +
                 Scalar(1.0) * _tmp107 * (-_tmp108 * _tmp129 - _tmp127 * _tmp128 + Scalar(1.0)) +
                 Scalar(1.0) * _tmp114 * (-_tmp128 * _tmp140 + _tmp139) +
                 Scalar(1.0) * _tmp120 * (_tmp117 * _tmp129 - _tmp128 * _tmp138) +
                 Scalar(1.0) * _tmp131 * (-_tmp128 * _tmp132 + _tmp133) +
                 Scalar(1.0) * _tmp134 * (_tmp129 * _tmp135 - _tmp129)));
  const Scalar _tmp142 = Scalar(1.4083112389913199) * _tmp121;
  const Scalar _tmp143 = _tmp129 * _tmp134;
  const Scalar _tmp144 = _tmp107 * _tmp109;
  const Scalar _tmp145 = _tmp100 * _tmp106 * _tmp75 - _tmp112 * _tmp114 +
                         _tmp119 * _tmp120 * _tmp75 - _tmp144 * _tmp76;
  const Scalar _tmp146 = Scalar(1.0) / (_tmp145);
  const Scalar _tmp147 =
      std::asinh(_tmp146 * (_tmp106 * _tmp137 * _tmp80 + _tmp107 * _tmp127 * _tmp80 +
                            _tmp114 * _tmp140 * _tmp80 + _tmp120 * _tmp138 * _tmp80 +
                            _tmp131 * _tmp132 * _tmp80 - _tmp135 * _tmp143));
  const Scalar _tmp148 = Scalar(1.4083112389913199) * _tmp145;
  const Scalar _tmp149 = _tmp106 * _tmp99 + _tmp110 * _tmp114 + _tmp118 * _tmp120 + _tmp144;
  const Scalar _tmp150 = Scalar(1.0) / (_tmp149);
  const Scalar _tmp151 =
      std::asinh(_tmp150 * (-_tmp106 * _tmp125 * _tmp98 + _tmp107 * _tmp126 - _tmp114 * _tmp139 -
                            _tmp117 * _tmp120 * _tmp125 - _tmp131 * _tmp133 + _tmp143));
  const Scalar _tmp152 = Scalar(1.4083112389913199) * _tmp149;

  // Output terms (1)
  Eigen::Matrix<Scalar, 4, 1> _res;

  _res(0, 0) =
      -_tmp30 * (Scalar(34.083374946563197) * _tmp0 + std::cosh(Scalar(1.0) * _tmp29) -
                 std::cosh(Scalar(0.71007031138673404) * _tmp0 *
                           (-_tmp29 * _tmp30 -
                            Scalar(125.0) *
                                std::sqrt(Scalar(
                                    Scalar(0.77439999999999998) *
                                        std::pow(Scalar(1 - Scalar(0.0090909090909090905) * _tmp28),
                                                 Scalar(2)) +
                                    std::pow(Scalar(-Scalar(0.0080000000000000002) * _tmp18 - 1),
                                             Scalar(2))))))) +
      _tmp35 + position_vector(2, 0);
  _res(1, 0) =
      -_tmp142 *
          (Scalar(34.083374946563197) * _tmp122 + std::cosh(Scalar(1.0) * _tmp141) -
           std::cosh(Scalar(0.71007031138673404) * _tmp122 *
                     (-_tmp141 * _tmp142 -
                      Scalar(125.0) *
                          std::sqrt(Scalar(
                              Scalar(0.77439999999999998) *
                                  std::pow(Scalar(-Scalar(0.0090909090909090905) * _tmp38 - 1),
                                           Scalar(2)) +
                              std::pow(Scalar(-Scalar(0.0080000000000000002) * _tmp42 - 1),
                                       Scalar(2))))))) +
      _tmp59 + position_vector(2, 0);
  _res(2, 0) =
      -_tmp148 *
          (Scalar(34.083374946563197) * _tmp146 + std::cosh(Scalar(1.0) * _tmp147) -
           std::cosh(
               Scalar(0.71007031138673404) * _tmp146 *
               (-_tmp147 * _tmp148 -
                Scalar(125.0) *
                    std::sqrt(Scalar(
                        std::pow(Scalar(1 - Scalar(0.0080000000000000002) * _tmp68), Scalar(2)) +
                        Scalar(0.77439999999999998) *
                            std::pow(Scalar(-Scalar(0.0090909090909090905) * _tmp65 - 1),
                                     Scalar(2))))))) +
      _tmp63 + position_vector(2, 0);
  _res(3, 0) =
      -_tmp152 *
          (Scalar(34.083374946563197) * _tmp150 + std::cosh(Scalar(1.0) * _tmp151) -
           std::cosh(
               Scalar(0.71007031138673404) * _tmp150 *
               (-_tmp151 * _tmp152 -
                Scalar(125.0) * std::sqrt(Scalar(
                                    Scalar(0.77439999999999998) *
                                        std::pow(Scalar(1 - Scalar(0.0090909090909090905) * _tmp48),
                                                 Scalar(2)) +
                                    std::pow(Scalar(1 - Scalar(0.0080000000000000002) * _tmp52),
                                             Scalar(2))))))) +
      _tmp62 + position_vector(2, 0);

  return _res;
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym

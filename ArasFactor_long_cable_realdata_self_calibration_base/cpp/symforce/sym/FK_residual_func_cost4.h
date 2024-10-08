// -----------------------------------------------------------------------------
// This file was autogenerated by symforce from template:
//     function/FUNCTION.h.jinja
// Do NOT modify by hand.
// -----------------------------------------------------------------------------

#pragma once

#include <Eigen/Dense>

#include <sym/pose3.h>
#include <sym/rot3.h>

namespace sym {

/**
 * This function was autogenerated from a symbolic function. Do not modify by hand.
 *
 * Symbolic function: FK_residual_func_cost4
 *
 * Args:
 *     fh1: Scalar
 *     fv1: Scalar
 *     DeltaRot: Rot3
 *     TransformationMatrix: Pose3
 *     encoder: Matrix41
 *     p_a: Matrix31
 *     p_b: Matrix31
 *     p_c: Matrix31
 *     p_d: Matrix31
 *     epsilon: Scalar
 *
 * Outputs:
 *     res: Matrix41
 */
template <typename Scalar>
Eigen::Matrix<Scalar, 4, 1> FkResidualFuncCost4(
    const Scalar fh1, const Scalar fv1, const sym::Rot3<Scalar>& DeltaRot,
    const sym::Pose3<Scalar>& TransformationMatrix, const Eigen::Matrix<Scalar, 4, 1>& encoder,
    const Eigen::Matrix<Scalar, 3, 1>& p_a, const Eigen::Matrix<Scalar, 3, 1>& p_b,
    const Eigen::Matrix<Scalar, 3, 1>& p_c, const Eigen::Matrix<Scalar, 3, 1>& p_d,
    const Scalar epsilon) {
  // Total ops: 460

  // Unused inputs
  (void)encoder;
  (void)epsilon;

  // Input arrays
  const Eigen::Matrix<Scalar, 4, 1>& _DeltaRot = DeltaRot.Data();
  const Eigen::Matrix<Scalar, 7, 1>& _TransformationMatrix = TransformationMatrix.Data();

  // Intermediate terms (156)
  const Scalar _tmp0 = Scalar(1.0) / (fh1);
  const Scalar _tmp1 = Scalar(0.71007031138673404) * _tmp0;
  const Scalar _tmp2 =
      _DeltaRot[0] * _TransformationMatrix[3] - _DeltaRot[1] * _TransformationMatrix[2] +
      _DeltaRot[2] * _TransformationMatrix[1] + _DeltaRot[3] * _TransformationMatrix[0];
  const Scalar _tmp3 = -2 * std::pow(_tmp2, Scalar(2));
  const Scalar _tmp4 =
      -_DeltaRot[0] * _TransformationMatrix[1] + _DeltaRot[1] * _TransformationMatrix[0] +
      _DeltaRot[2] * _TransformationMatrix[3] + _DeltaRot[3] * _TransformationMatrix[2];
  const Scalar _tmp5 = 1 - 2 * std::pow(_tmp4, Scalar(2));
  const Scalar _tmp6 = Scalar(0.20999999999999999) * _tmp3 + Scalar(0.20999999999999999) * _tmp5;
  const Scalar _tmp7 = -_tmp6;
  const Scalar _tmp8 =
      _DeltaRot[0] * _TransformationMatrix[2] + _DeltaRot[1] * _TransformationMatrix[3] -
      _DeltaRot[2] * _TransformationMatrix[0] + _DeltaRot[3] * _TransformationMatrix[1];
  const Scalar _tmp9 = 2 * _tmp2;
  const Scalar _tmp10 = _tmp8 * _tmp9;
  const Scalar _tmp11 =
      -2 * _DeltaRot[0] * _TransformationMatrix[0] - 2 * _DeltaRot[1] * _TransformationMatrix[1] -
      2 * _DeltaRot[2] * _TransformationMatrix[2] + 2 * _DeltaRot[3] * _TransformationMatrix[3];
  const Scalar _tmp12 = _tmp11 * _tmp4;
  const Scalar _tmp13 = Scalar(0.20999999999999999) * _tmp10 + Scalar(0.20999999999999999) * _tmp12;
  const Scalar _tmp14 = 2 * _tmp4 * _tmp8;
  const Scalar _tmp15 = _tmp11 * _tmp2;
  const Scalar _tmp16 = _tmp14 - _tmp15;
  const Scalar _tmp17 = -Scalar(0.010999999999999999) * _tmp16;
  const Scalar _tmp18 = -_tmp13 + _tmp17;
  const Scalar _tmp19 = _tmp18 + _tmp7;
  const Scalar _tmp20 = _TransformationMatrix[5] + _tmp19;
  const Scalar _tmp21 = -2 * std::pow(_tmp8, Scalar(2));
  const Scalar _tmp22 = Scalar(0.20999999999999999) * _tmp21 + Scalar(0.20999999999999999) * _tmp5;
  const Scalar _tmp23 = -_tmp22;
  const Scalar _tmp24 = _tmp4 * _tmp9;
  const Scalar _tmp25 = _tmp11 * _tmp8;
  const Scalar _tmp26 = _tmp24 + _tmp25;
  const Scalar _tmp27 = -Scalar(0.010999999999999999) * _tmp26;
  const Scalar _tmp28 = Scalar(0.20999999999999999) * _tmp10 - Scalar(0.20999999999999999) * _tmp12;
  const Scalar _tmp29 = _tmp27 - _tmp28;
  const Scalar _tmp30 = _tmp23 + _tmp29;
  const Scalar _tmp31 = _TransformationMatrix[4] + _tmp30;
  const Scalar _tmp32 = Scalar(1.4083112389913199) * fh1;
  const Scalar _tmp33 =
      std::cosh(_tmp1 * (-_tmp32 * std::asinh(_tmp0 * fv1) -
                         std::sqrt(Scalar(std::pow(Scalar(-_tmp20 + p_a(1, 0)), Scalar(2)) +
                                          std::pow(Scalar(-_tmp31 + p_a(0, 0)), Scalar(2))))));
  const Scalar _tmp34 = Scalar(333.54000000000002) - fv1;
  const Scalar _tmp35 = Scalar(0.20999999999999999) * _tmp24 - Scalar(0.20999999999999999) * _tmp25;
  const Scalar _tmp36 = -_tmp35;
  const Scalar _tmp37 = -Scalar(0.010999999999999999) * _tmp21 -
                        Scalar(0.010999999999999999) * _tmp3 + Scalar(-0.010999999999999999);
  const Scalar _tmp38 = Scalar(0.20999999999999999) * _tmp14 + Scalar(0.20999999999999999) * _tmp15;
  const Scalar _tmp39 = _tmp37 + _tmp38;
  const Scalar _tmp40 = _tmp36 + _tmp39;
  const Scalar _tmp41 = _tmp18 + _tmp6;
  const Scalar _tmp42 = _TransformationMatrix[5] + _tmp41;
  const Scalar _tmp43 = _tmp42 - p_d(1, 0);
  const Scalar _tmp44 = _tmp27 + _tmp28;
  const Scalar _tmp45 = _tmp23 + _tmp44;
  const Scalar _tmp46 = _TransformationMatrix[4] + _tmp45;
  const Scalar _tmp47 = _tmp46 - p_d(0, 0);
  const Scalar _tmp48 = std::pow(Scalar(std::pow(_tmp43, Scalar(2)) + std::pow(_tmp47, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp49 = _tmp47 * _tmp48;
  const Scalar _tmp50 = _tmp37 - _tmp38;
  const Scalar _tmp51 = _tmp35 + _tmp50;
  const Scalar _tmp52 = _tmp35 + _tmp39;
  const Scalar _tmp53 = _tmp22 + _tmp44;
  const Scalar _tmp54 = _TransformationMatrix[4] + _tmp53;
  const Scalar _tmp55 = _tmp54 - p_c(0, 0);
  const Scalar _tmp56 = _tmp13 + _tmp17;
  const Scalar _tmp57 = _tmp56 + _tmp6;
  const Scalar _tmp58 = _TransformationMatrix[5] + _tmp57;
  const Scalar _tmp59 = _tmp58 - p_c(1, 0);
  const Scalar _tmp60 = std::pow(Scalar(std::pow(_tmp55, Scalar(2)) + std::pow(_tmp59, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp61 = _tmp55 * _tmp60;
  const Scalar _tmp62 = _tmp51 * _tmp61 - _tmp52 * _tmp61;
  const Scalar _tmp63 = _tmp43 * _tmp48;
  const Scalar _tmp64 = _tmp56 + _tmp7;
  const Scalar _tmp65 = _TransformationMatrix[5] + _tmp64;
  const Scalar _tmp66 = _tmp65 - p_b(1, 0);
  const Scalar _tmp67 = _tmp22 + _tmp29;
  const Scalar _tmp68 = _TransformationMatrix[4] + _tmp67;
  const Scalar _tmp69 = _tmp68 - p_b(0, 0);
  const Scalar _tmp70 = Scalar(1.0) / (_tmp69);
  const Scalar _tmp71 = _tmp66 * _tmp70;
  const Scalar _tmp72 = _tmp49 * _tmp71 - _tmp63;
  const Scalar _tmp73 = _tmp59 * _tmp60;
  const Scalar _tmp74 = Scalar(1.0) / (_tmp61 * _tmp71 - _tmp73);
  const Scalar _tmp75 = _tmp72 * _tmp74;
  const Scalar _tmp76 = _tmp51 * _tmp71;
  const Scalar _tmp77 = _tmp74 * (_tmp52 * _tmp73 - _tmp61 * _tmp76);
  const Scalar _tmp78 = _tmp40 * _tmp63 - _tmp49 * _tmp76 - _tmp72 * _tmp77;
  const Scalar _tmp79 = Scalar(1.0) * _tmp64;
  const Scalar _tmp80 = -_tmp79;
  const Scalar _tmp81 = Scalar(1.0) / (_tmp57 + _tmp80);
  const Scalar _tmp82 = Scalar(1.0) * _tmp67;
  const Scalar _tmp83 = -_tmp53 + _tmp82;
  const Scalar _tmp84 = _tmp81 * _tmp83;
  const Scalar _tmp85 = -_tmp40 * _tmp49 + _tmp49 * _tmp51 - _tmp62 * _tmp75 - _tmp78 * _tmp84;
  const Scalar _tmp86 = Scalar(1.0) / (_tmp85);
  const Scalar _tmp87 = _tmp79 * _tmp84 + _tmp82;
  const Scalar _tmp88 = 0;
  const Scalar _tmp89 = _tmp61 * _tmp75;
  const Scalar _tmp90 =
      std::sqrt(Scalar(std::pow(_tmp66, Scalar(2)) + std::pow(_tmp69, Scalar(2))));
  const Scalar _tmp91 = _tmp70 * _tmp90;
  const Scalar _tmp92 = Scalar(1.0) * _tmp74;
  const Scalar _tmp93 = Scalar(1.0) * _tmp81;
  const Scalar _tmp94 = -_tmp62 * _tmp92 + _tmp77 * _tmp83 * _tmp93;
  const Scalar _tmp95 = Scalar(1.0) / (_tmp90);
  const Scalar _tmp96 = _tmp91 * (-_tmp64 * _tmp69 * _tmp95 + _tmp66 * _tmp67 * _tmp95);
  const Scalar _tmp97 = -_tmp53 * _tmp73 + _tmp57 * _tmp61 + _tmp61 * _tmp96;
  const Scalar _tmp98 = _tmp41 * _tmp49 - _tmp45 * _tmp63 + _tmp49 * _tmp96 - _tmp75 * _tmp97;
  const Scalar _tmp99 = _tmp86 * _tmp98;
  const Scalar _tmp100 = Scalar(1.0) / (_tmp98);
  const Scalar _tmp101 = _tmp100 * _tmp85;
  const Scalar _tmp102 = _tmp101 * (-_tmp92 * _tmp97 - _tmp94 * _tmp99);
  const Scalar _tmp103 = _tmp86 * (_tmp102 + _tmp94);
  const Scalar _tmp104 = -_tmp103 * _tmp72 + Scalar(1.0);
  const Scalar _tmp105 = _tmp61 * _tmp74;
  const Scalar _tmp106 = _tmp31 - p_a(0, 0);
  const Scalar _tmp107 = _tmp20 - p_a(1, 0);
  const Scalar _tmp108 =
      std::pow(Scalar(std::pow(_tmp106, Scalar(2)) + std::pow(_tmp107, Scalar(2))),
               Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp109 = _tmp107 * _tmp108;
  const Scalar _tmp110 = _tmp109 * fh1;
  const Scalar _tmp111 = _tmp71 * _tmp74;
  const Scalar _tmp112 = _tmp71 * _tmp77 + _tmp76;
  const Scalar _tmp113 = _tmp111 * _tmp62 - _tmp112 * _tmp84 - _tmp51;
  const Scalar _tmp114 = _tmp101 * (_tmp111 * _tmp97 - _tmp113 * _tmp99 - _tmp96);
  const Scalar _tmp115 = _tmp86 * (_tmp113 + _tmp114);
  const Scalar _tmp116 = -_tmp115 * _tmp72 - _tmp71;
  const Scalar _tmp117 = _tmp106 * _tmp108;
  const Scalar _tmp118 = _tmp117 * fh1;
  const Scalar _tmp119 = Scalar(1.0) * _tmp100;
  const Scalar _tmp120 = fh1 * (_tmp109 * _tmp30 - _tmp117 * _tmp19);
  const Scalar _tmp121 = -_tmp110 * _tmp91 * (_tmp103 * _tmp49 + _tmp104 * _tmp105) -
                         _tmp118 * _tmp91 * (_tmp105 * _tmp116 + _tmp115 * _tmp49 + Scalar(1.0)) -
                         _tmp120 * _tmp91 * (_tmp119 * _tmp49 - _tmp119 * _tmp89) -
                         _tmp34 * _tmp91 * (_tmp49 * _tmp88 - _tmp88 * _tmp89);
  const Scalar _tmp122 = Scalar(1.0) / (_tmp121);
  const Scalar _tmp123 = Scalar(0.71007031138673404) * _tmp122;
  const Scalar _tmp124 = fh1 * (_tmp36 + _tmp50);
  const Scalar _tmp125 = -_tmp109 * _tmp124 - Scalar(40.024799999999999) * _tmp16 - _tmp19 * fv1;
  const Scalar _tmp126 = _tmp41 + _tmp80;
  const Scalar _tmp127 = _tmp126 * _tmp84;
  const Scalar _tmp128 = Scalar(1.0) / (-_tmp127 - _tmp45 + _tmp82);
  const Scalar _tmp129 = Scalar(1.0) * _tmp128;
  const Scalar _tmp130 = _tmp129 * _tmp84;
  const Scalar _tmp131 = _tmp127 * _tmp129 + Scalar(1.0);
  const Scalar _tmp132 = _tmp101 * _tmp129;
  const Scalar _tmp133 = -_tmp119 * _tmp78 + _tmp126 * _tmp132;
  const Scalar _tmp134 = _tmp126 * _tmp128;
  const Scalar _tmp135 = _tmp112 + _tmp114 * _tmp134 - _tmp115 * _tmp78;
  const Scalar _tmp136 = _tmp102 * _tmp134 - _tmp103 * _tmp78 - Scalar(1.0) * _tmp77;
  const Scalar _tmp137 = _tmp126 * _tmp81;
  const Scalar _tmp138 = _tmp117 * _tmp124 + Scalar(40.024799999999999) * _tmp26 + _tmp30 * fv1;
  const Scalar _tmp139 = _tmp128 * _tmp87;
  const Scalar _tmp140 = -_tmp126 * _tmp139 - _tmp78 * _tmp88 + _tmp80;
  const Scalar _tmp141 = Scalar(1.4083112389913199) * _tmp121;
  const Scalar _tmp142 = std::cosh(
      _tmp123 *
      (-_tmp141 *
           std::asinh(_tmp122 * (Scalar(1.0) * _tmp110 * (_tmp102 * _tmp129 - _tmp136 * _tmp93) +
                                 Scalar(1.0) * _tmp118 * (_tmp114 * _tmp129 - _tmp135 * _tmp93) +
                                 Scalar(1.0) * _tmp120 * (_tmp132 - _tmp133 * _tmp93) +
                                 Scalar(1.0) * _tmp125 * (_tmp130 - _tmp131 * _tmp93) +
                                 Scalar(1.0) * _tmp138 * (_tmp129 * _tmp137 - _tmp129) +
                                 Scalar(1.0) * _tmp34 *
                                     (-_tmp129 * _tmp87 - _tmp140 * _tmp93 + Scalar(1.0)))) -
       std::sqrt(Scalar(std::pow(Scalar(-_tmp65 + p_b(1, 0)), Scalar(2)) +
                        std::pow(Scalar(-_tmp68 + p_b(0, 0)), Scalar(2))))));
  const Scalar _tmp143 = _tmp119 * _tmp120;
  const Scalar _tmp144 = _tmp34 * _tmp88;
  const Scalar _tmp145 =
      _tmp104 * _tmp110 * _tmp74 + _tmp116 * _tmp118 * _tmp74 - _tmp143 * _tmp75 - _tmp144 * _tmp75;
  const Scalar _tmp146 = Scalar(1.0) / (_tmp145);
  const Scalar _tmp147 = Scalar(0.71007031138673404) * _tmp146;
  const Scalar _tmp148 = _tmp129 * _tmp138;
  const Scalar _tmp149 = Scalar(1.4083112389913199) * _tmp145;
  const Scalar _tmp150 = std::cosh(
      _tmp147 *
      (-_tmp149 * std::asinh(_tmp146 * (_tmp110 * _tmp136 * _tmp81 + _tmp118 * _tmp135 * _tmp81 +
                                        _tmp120 * _tmp133 * _tmp81 + _tmp125 * _tmp131 * _tmp81 -
                                        _tmp137 * _tmp148 + _tmp140 * _tmp34 * _tmp81)) -
       std::sqrt(Scalar(std::pow(Scalar(-_tmp54 + p_c(0, 0)), Scalar(2)) +
                        std::pow(Scalar(-_tmp58 + p_c(1, 0)), Scalar(2))))));
  const Scalar _tmp151 = _tmp103 * _tmp110 + _tmp115 * _tmp118 + _tmp143 + _tmp144;
  const Scalar _tmp152 = Scalar(1.0) / (_tmp151);
  const Scalar _tmp153 = Scalar(0.71007031138673404) * _tmp152;
  const Scalar _tmp154 = Scalar(1.4083112389913199) * _tmp151;
  const Scalar _tmp155 = std::cosh(
      _tmp153 * (-_tmp154 * std::asinh(_tmp152 * (-_tmp102 * _tmp110 * _tmp128 -
                                                  _tmp114 * _tmp118 * _tmp128 - _tmp120 * _tmp132 -
                                                  _tmp125 * _tmp130 + _tmp139 * _tmp34 + _tmp148)) -
                 std::sqrt(Scalar(std::pow(Scalar(-_tmp42 + p_d(1, 0)), Scalar(2)) +
                                  std::pow(Scalar(-_tmp46 + p_d(0, 0)), Scalar(2))))));

  // Output terms (1)
  Eigen::Matrix<Scalar, 4, 1> _res;

  _res(0, 0) = -_tmp1 * _tmp32 * p_a(2, 0) + p_a(2, 0);
  _res(1, 0) = -_tmp123 * _tmp141 * p_b(2, 0) + p_b(2, 0);
  _res(2, 0) = -_tmp147 * _tmp149 * p_c(2, 0) + p_c(2, 0);
  _res(3, 0) = -_tmp153 * _tmp154 * p_d(2, 0) + p_d(2, 0);

  return _res;
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym

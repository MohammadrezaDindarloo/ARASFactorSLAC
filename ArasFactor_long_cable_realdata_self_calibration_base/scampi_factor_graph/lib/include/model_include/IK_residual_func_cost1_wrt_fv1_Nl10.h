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
 * Symbolic function: IK_residual_func_cost1_wrt_fv1_Nl10
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
Eigen::Matrix<Scalar, 4, 1> IkResidualFuncCost1WrtFv1Nl10(
    const Scalar fh1, const Scalar fv1, const sym::Rot3<Scalar>& DeltaRot,
    const Eigen::Matrix<Scalar, 3, 1>& position_vector, const Eigen::Matrix<Scalar, 4, 1>& encoder,
    const Eigen::Matrix<Scalar, 3, 1>& p_a, const Eigen::Matrix<Scalar, 3, 1>& p_b,
    const Eigen::Matrix<Scalar, 3, 1>& p_c, const Eigen::Matrix<Scalar, 3, 1>& p_d,
    const sym::Rot3<Scalar>& Rot_init, const Scalar epsilon) {
  // Total ops: 586

  // Unused inputs
  (void)encoder;
  (void)epsilon;

  // Input arrays
  const Eigen::Matrix<Scalar, 4, 1>& _DeltaRot = DeltaRot.Data();
  const Eigen::Matrix<Scalar, 4, 1>& _Rot_init = Rot_init.Data();

  // Intermediate terms (195)
  const Scalar _tmp0 = Scalar(1.0) / (fh1);
  const Scalar _tmp1 = std::asinh(_tmp0 * fv1);
  const Scalar _tmp2 = Scalar(1.0) * _tmp0 /
                       std::sqrt(Scalar(1 + std::pow(fv1, Scalar(2)) / std::pow(fh1, Scalar(2))));
  const Scalar _tmp3 = _DeltaRot[0] * _Rot_init[3] - _DeltaRot[1] * _Rot_init[2] +
                       _DeltaRot[2] * _Rot_init[1] + _DeltaRot[3] * _Rot_init[0];
  const Scalar _tmp4 = _DeltaRot[0] * _Rot_init[2] + _DeltaRot[1] * _Rot_init[3] -
                       _DeltaRot[2] * _Rot_init[0] + _DeltaRot[3] * _Rot_init[1];
  const Scalar _tmp5 = 2 * _tmp4;
  const Scalar _tmp6 = _tmp3 * _tmp5;
  const Scalar _tmp7 = -_DeltaRot[0] * _Rot_init[1] + _DeltaRot[1] * _Rot_init[0] +
                       _DeltaRot[2] * _Rot_init[3] + _DeltaRot[3] * _Rot_init[2];
  const Scalar _tmp8 = -2 * _DeltaRot[0] * _Rot_init[0] - 2 * _DeltaRot[1] * _Rot_init[1] -
                       2 * _DeltaRot[2] * _Rot_init[2] + 2 * _DeltaRot[3] * _Rot_init[3];
  const Scalar _tmp9 = _tmp7 * _tmp8;
  const Scalar _tmp10 = Scalar(0.20999999999999999) * _tmp6 + Scalar(0.20999999999999999) * _tmp9;
  const Scalar _tmp11 = -2 * std::pow(_tmp3, Scalar(2));
  const Scalar _tmp12 = -2 * std::pow(_tmp7, Scalar(2));
  const Scalar _tmp13 = Scalar(0.20999999999999999) * _tmp11 +
                        Scalar(0.20999999999999999) * _tmp12 + Scalar(0.20999999999999999);
  const Scalar _tmp14 = _tmp5 * _tmp7;
  const Scalar _tmp15 = _tmp3 * _tmp8;
  const Scalar _tmp16 = _tmp14 - _tmp15;
  const Scalar _tmp17 = Scalar(0.010999999999999999) * _tmp16;
  const Scalar _tmp18 = -_tmp17;
  const Scalar _tmp19 = -_tmp13 + _tmp18;
  const Scalar _tmp20 = _tmp10 + _tmp19;
  const Scalar _tmp21 = _tmp20 + position_vector(1, 0);
  const Scalar _tmp22 = 1 - 2 * std::pow(_tmp4, Scalar(2));
  const Scalar _tmp23 = Scalar(0.20999999999999999) * _tmp12 + Scalar(0.20999999999999999) * _tmp22;
  const Scalar _tmp24 = 2 * _tmp3 * _tmp7;
  const Scalar _tmp25 = _tmp4 * _tmp8;
  const Scalar _tmp26 = _tmp24 + _tmp25;
  const Scalar _tmp27 = -Scalar(0.010999999999999999) * _tmp26;
  const Scalar _tmp28 = Scalar(0.20999999999999999) * _tmp6 - Scalar(0.20999999999999999) * _tmp9;
  const Scalar _tmp29 = _tmp27 - _tmp28;
  const Scalar _tmp30 = _tmp23 + _tmp29;
  const Scalar _tmp31 = _tmp30 + position_vector(0, 0);
  const Scalar _tmp32 = Scalar(1.4083112389913199) * fh1;
  const Scalar _tmp33 = Scalar(333.54000000000002) - fv1;
  const Scalar _tmp34 = Scalar(0.20999999999999999) * _tmp24 - Scalar(0.20999999999999999) * _tmp25;
  const Scalar _tmp35 =
      -Scalar(0.010999999999999999) * _tmp11 - Scalar(0.010999999999999999) * _tmp22;
  const Scalar _tmp36 = Scalar(0.20999999999999999) * _tmp14 + Scalar(0.20999999999999999) * _tmp15;
  const Scalar _tmp37 = _tmp35 + _tmp36;
  const Scalar _tmp38 = _tmp34 + _tmp37;
  const Scalar _tmp39 = _tmp27 + _tmp28;
  const Scalar _tmp40 = _tmp23 + _tmp39;
  const Scalar _tmp41 = _tmp40 + position_vector(0, 0);
  const Scalar _tmp42 = _tmp41 - p_c(0, 0);
  const Scalar _tmp43 = _tmp10 + _tmp13 + _tmp18;
  const Scalar _tmp44 = _tmp43 + position_vector(1, 0);
  const Scalar _tmp45 = _tmp44 - p_c(1, 0);
  const Scalar _tmp46 = std::pow(Scalar(std::pow(_tmp42, Scalar(2)) + std::pow(_tmp45, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp47 = _tmp42 * _tmp46;
  const Scalar _tmp48 = -_tmp34;
  const Scalar _tmp49 = _tmp35 - _tmp36;
  const Scalar _tmp50 = _tmp48 + _tmp49;
  const Scalar _tmp51 = -_tmp23;
  const Scalar _tmp52 = _tmp29 + _tmp51;
  const Scalar _tmp53 = _tmp52 + position_vector(0, 0);
  const Scalar _tmp54 = _tmp53 - p_a(0, 0);
  const Scalar _tmp55 = -_tmp10;
  const Scalar _tmp56 = _tmp19 + _tmp55;
  const Scalar _tmp57 = _tmp56 + position_vector(1, 0);
  const Scalar _tmp58 = _tmp57 - p_a(1, 0);
  const Scalar _tmp59 = std::pow(Scalar(std::pow(_tmp54, Scalar(2)) + std::pow(_tmp58, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp60 = _tmp54 * _tmp59;
  const Scalar _tmp61 = _tmp37 + _tmp48;
  const Scalar _tmp62 = _tmp60 * _tmp61;
  const Scalar _tmp63 = -_tmp50 * _tmp60 + _tmp62;
  const Scalar _tmp64 = _tmp39 + _tmp51;
  const Scalar _tmp65 = _tmp64 + position_vector(0, 0);
  const Scalar _tmp66 = _tmp65 - p_d(0, 0);
  const Scalar _tmp67 = Scalar(1.0) / (_tmp66);
  const Scalar _tmp68 = _tmp13 + _tmp55;
  const Scalar _tmp69 = _tmp18 + _tmp68;
  const Scalar _tmp70 = _tmp69 + position_vector(1, 0);
  const Scalar _tmp71 = _tmp70 - p_d(1, 0);
  const Scalar _tmp72 = _tmp67 * _tmp71;
  const Scalar _tmp73 = _tmp58 * _tmp59;
  const Scalar _tmp74 = Scalar(1.0) / (_tmp60 * _tmp72 - _tmp73);
  const Scalar _tmp75 = _tmp45 * _tmp46;
  const Scalar _tmp76 = _tmp47 * _tmp72 - _tmp75;
  const Scalar _tmp77 = _tmp74 * _tmp76;
  const Scalar _tmp78 = _tmp47 * _tmp61;
  const Scalar _tmp79 = _tmp50 * _tmp73 - _tmp62 * _tmp72;
  const Scalar _tmp80 = _tmp38 * _tmp75 - _tmp72 * _tmp78 - _tmp77 * _tmp79;
  const Scalar _tmp81 = Scalar(1.0) * _tmp69;
  const Scalar _tmp82 = -_tmp81;
  const Scalar _tmp83 = Scalar(1.0) / (_tmp56 + _tmp82);
  const Scalar _tmp84 = Scalar(1.0) * _tmp64;
  const Scalar _tmp85 = _tmp83 * (-_tmp52 + _tmp84);
  const Scalar _tmp86 = -_tmp38 * _tmp47 - _tmp63 * _tmp77 + _tmp78 - _tmp80 * _tmp85;
  const Scalar _tmp87 = Scalar(1.0) / (_tmp86);
  const Scalar _tmp88 = _tmp81 * _tmp85 + _tmp84;
  const Scalar _tmp89 = 0;
  const Scalar _tmp90 = _tmp87 * _tmp89;
  const Scalar _tmp91 = _tmp60 * _tmp77;
  const Scalar _tmp92 =
      std::sqrt(Scalar(std::pow(_tmp66, Scalar(2)) + std::pow(_tmp71, Scalar(2))));
  const Scalar _tmp93 = _tmp67 * _tmp92;
  const Scalar _tmp94 = _tmp93 * (_tmp47 * _tmp90 - _tmp90 * _tmp91);
  const Scalar _tmp95 = Scalar(1.0) * _tmp74;
  const Scalar _tmp96 = _tmp79 * _tmp95;
  const Scalar _tmp97 = -_tmp63 * _tmp95 + _tmp85 * _tmp96;
  const Scalar _tmp98 = Scalar(1.0) / (_tmp92);
  const Scalar _tmp99 = _tmp93 * (_tmp64 * _tmp71 * _tmp98 - _tmp66 * _tmp69 * _tmp98);
  const Scalar _tmp100 = -_tmp52 * _tmp73 + _tmp56 * _tmp60 + _tmp60 * _tmp99;
  const Scalar _tmp101 = -_tmp100 * _tmp77 - _tmp40 * _tmp75 + _tmp43 * _tmp47 + _tmp47 * _tmp99;
  const Scalar _tmp102 = _tmp101 * _tmp87;
  const Scalar _tmp103 = Scalar(1.0) / (_tmp101);
  const Scalar _tmp104 = _tmp103 * _tmp86;
  const Scalar _tmp105 = _tmp104 * (-_tmp100 * _tmp95 - _tmp102 * _tmp97);
  const Scalar _tmp106 = _tmp105 + _tmp97;
  const Scalar _tmp107 = _tmp47 * _tmp87;
  const Scalar _tmp108 = _tmp76 * _tmp87;
  const Scalar _tmp109 = -_tmp106 * _tmp108 + Scalar(1.0);
  const Scalar _tmp110 = _tmp60 * _tmp74;
  const Scalar _tmp111 = _tmp31 - p_b(0, 0);
  const Scalar _tmp112 = _tmp21 - p_b(1, 0);
  const Scalar _tmp113 =
      std::pow(Scalar(std::pow(_tmp111, Scalar(2)) + std::pow(_tmp112, Scalar(2))),
               Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp114 = _tmp112 * _tmp113;
  const Scalar _tmp115 = _tmp114 * fh1;
  const Scalar _tmp116 = _tmp72 * _tmp74;
  const Scalar _tmp117 = _tmp116 * _tmp79 + _tmp61 * _tmp72;
  const Scalar _tmp118 = _tmp116 * _tmp63 - _tmp117 * _tmp85 - _tmp61;
  const Scalar _tmp119 = _tmp104 * (_tmp100 * _tmp116 - _tmp102 * _tmp118 - _tmp99);
  const Scalar _tmp120 = _tmp118 + _tmp119;
  const Scalar _tmp121 = -_tmp108 * _tmp120 - _tmp72;
  const Scalar _tmp122 = _tmp111 * _tmp113;
  const Scalar _tmp123 = _tmp122 * fh1;
  const Scalar _tmp124 = Scalar(1.0) * _tmp103;
  const Scalar _tmp125 = fh1 * (_tmp114 * _tmp30 - _tmp122 * _tmp20);
  const Scalar _tmp126 = -_tmp115 * _tmp93 * (_tmp106 * _tmp107 + _tmp109 * _tmp110) -
                         _tmp123 * _tmp93 * (_tmp107 * _tmp120 + _tmp110 * _tmp121 + Scalar(1.0)) -
                         _tmp125 * _tmp93 * (_tmp124 * _tmp47 - _tmp124 * _tmp91) - _tmp33 * _tmp94;
  const Scalar _tmp127 = std::pow(_tmp126, Scalar(-2));
  const Scalar _tmp128 = _tmp127 * _tmp94;
  const Scalar _tmp129 = Scalar(0.71007031138673404) * _tmp128;
  const Scalar _tmp130 = Scalar(1.0) / (_tmp126);
  const Scalar _tmp131 = fh1 * (_tmp34 + _tmp49);
  const Scalar _tmp132 = _tmp122 * _tmp131 + Scalar(40.024799999999999) * _tmp26 + _tmp30 * fv1;
  const Scalar _tmp133 = _tmp43 + _tmp82;
  const Scalar _tmp134 = _tmp133 * _tmp85;
  const Scalar _tmp135 = Scalar(1.0) / (-_tmp134 - _tmp40 + _tmp84);
  const Scalar _tmp136 = Scalar(1.0) * _tmp135;
  const Scalar _tmp137 = _tmp133 * _tmp83;
  const Scalar _tmp138 = Scalar(1.0) * _tmp136 * _tmp137 - Scalar(1.0) * _tmp136;
  const Scalar _tmp139 = -_tmp114 * _tmp131 - Scalar(40.024799999999999) * _tmp16 - _tmp20 * fv1;
  const Scalar _tmp140 = _tmp136 * _tmp85;
  const Scalar _tmp141 = _tmp83 * (_tmp134 * _tmp136 + Scalar(1.0));
  const Scalar _tmp142 = Scalar(1.0) * _tmp140 - Scalar(1.0) * _tmp141;
  const Scalar _tmp143 = _tmp104 * _tmp136;
  const Scalar _tmp144 = -_tmp124 * _tmp80 + _tmp133 * _tmp143;
  const Scalar _tmp145 = Scalar(1.0) * _tmp83;
  const Scalar _tmp146 = _tmp135 * _tmp88;
  const Scalar _tmp147 = _tmp80 * _tmp87;
  const Scalar _tmp148 = _tmp83 * (-_tmp133 * _tmp146 - _tmp147 * _tmp89 + _tmp82);
  const Scalar _tmp149 = -Scalar(1.0) * _tmp136 * _tmp88 - Scalar(1.0) * _tmp148 + Scalar(1.0);
  const Scalar _tmp150 = _tmp133 * _tmp135;
  const Scalar _tmp151 = _tmp117 + _tmp119 * _tmp150 - _tmp120 * _tmp147;
  const Scalar _tmp152 = _tmp105 * _tmp150 - _tmp106 * _tmp147 - _tmp96;
  const Scalar _tmp153 = Scalar(1.0) * _tmp115 * (_tmp105 * _tmp136 - _tmp145 * _tmp152) +
                         Scalar(1.0) * _tmp123 * (_tmp119 * _tmp136 - _tmp145 * _tmp151) +
                         Scalar(1.0) * _tmp125 * (_tmp143 - _tmp144 * _tmp145) + _tmp132 * _tmp138 +
                         _tmp139 * _tmp142 + _tmp149 * _tmp33;
  const Scalar _tmp154 = std::asinh(_tmp130 * _tmp153);
  const Scalar _tmp155 = Scalar(1.0) * _tmp154;
  const Scalar _tmp156 = _tmp17 + _tmp68;
  const Scalar _tmp157 =
      (-_tmp128 * _tmp153 + _tmp130 * (_tmp138 * _tmp30 + _tmp142 * _tmp156 - _tmp149)) /
      std::sqrt(Scalar(_tmp127 * std::pow(_tmp153, Scalar(2)) + 1));
  const Scalar _tmp158 = Scalar(1.4083112389913199) * _tmp126;
  const Scalar _tmp159 =
      -_tmp154 * _tmp158 - std::sqrt(Scalar(std::pow(Scalar(-_tmp65 + p_d(0, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp70 + p_d(1, 0)), Scalar(2))));
  const Scalar _tmp160 = Scalar(1.4083112389913199) * _tmp94;
  const Scalar _tmp161 = Scalar(0.71007031138673404) * _tmp130;
  const Scalar _tmp162 = _tmp159 * _tmp161;
  const Scalar _tmp163 = _tmp124 * _tmp125;
  const Scalar _tmp164 = _tmp33 * _tmp90;
  const Scalar _tmp165 =
      _tmp109 * _tmp115 * _tmp74 + _tmp121 * _tmp123 * _tmp74 - _tmp163 * _tmp77 - _tmp164 * _tmp77;
  const Scalar _tmp166 = Scalar(1.0) / (_tmp165);
  const Scalar _tmp167 = _tmp132 * _tmp136;
  const Scalar _tmp168 = _tmp115 * _tmp152 * _tmp83 + _tmp123 * _tmp151 * _tmp83 +
                         _tmp125 * _tmp144 * _tmp83 - _tmp137 * _tmp167 + _tmp139 * _tmp141 +
                         _tmp148 * _tmp33;
  const Scalar _tmp169 = std::asinh(_tmp166 * _tmp168);
  const Scalar _tmp170 = Scalar(1.0) * _tmp169;
  const Scalar _tmp171 = Scalar(0.71007031138673404) * _tmp166;
  const Scalar _tmp172 = Scalar(1.4083112389913199) * _tmp165;
  const Scalar _tmp173 =
      -_tmp169 * _tmp172 - std::sqrt(Scalar(std::pow(Scalar(-_tmp53 + p_a(0, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp57 + p_a(1, 0)), Scalar(2))));
  const Scalar _tmp174 = _tmp171 * _tmp173;
  const Scalar _tmp175 = Scalar(1.4083112389913199) * _tmp90;
  const Scalar _tmp176 = _tmp175 * _tmp77;
  const Scalar _tmp177 = std::pow(_tmp165, Scalar(-2));
  const Scalar _tmp178 = _tmp177 * _tmp77 * _tmp90;
  const Scalar _tmp179 = _tmp136 * _tmp30;
  const Scalar _tmp180 =
      (_tmp166 * (-_tmp137 * _tmp179 + _tmp141 * _tmp156 - _tmp148) - _tmp168 * _tmp178) /
      std::sqrt(Scalar(std::pow(_tmp168, Scalar(2)) * _tmp177 + 1));
  const Scalar _tmp181 = Scalar(0.71007031138673404) * _tmp178;
  const Scalar _tmp182 =
      _tmp106 * _tmp115 * _tmp87 + _tmp120 * _tmp123 * _tmp87 + _tmp163 + _tmp164;
  const Scalar _tmp183 = std::pow(_tmp182, Scalar(-2));
  const Scalar _tmp184 = _tmp183 * _tmp90;
  const Scalar _tmp185 = Scalar(0.71007031138673404) * _tmp184;
  const Scalar _tmp186 = -_tmp105 * _tmp115 * _tmp135 - _tmp119 * _tmp123 * _tmp135 -
                         _tmp125 * _tmp143 - _tmp139 * _tmp140 + _tmp146 * _tmp33 + _tmp167;
  const Scalar _tmp187 = Scalar(1.0) / (_tmp182);
  const Scalar _tmp188 = std::asinh(_tmp186 * _tmp187);
  const Scalar _tmp189 = Scalar(1.4083112389913199) * _tmp182;
  const Scalar _tmp190 =
      -_tmp188 * _tmp189 - std::sqrt(Scalar(std::pow(Scalar(-_tmp41 + p_c(0, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp44 + p_c(1, 0)), Scalar(2))));
  const Scalar _tmp191 = (_tmp184 * _tmp186 + _tmp187 * (-_tmp140 * _tmp156 - _tmp146 + _tmp179)) /
                         std::sqrt(Scalar(_tmp183 * std::pow(_tmp186, Scalar(2)) + 1));
  const Scalar _tmp192 = Scalar(0.71007031138673404) * _tmp187;
  const Scalar _tmp193 = _tmp190 * _tmp192;
  const Scalar _tmp194 = Scalar(1.0) * _tmp188;

  // Output terms (1)
  Eigen::Matrix<Scalar, 4, 1> _res;

  _res(0, 0) =
      -_tmp32 *
      (_tmp2 * std::sinh(Scalar(1.0) * _tmp1) +
       _tmp2 * std::sinh(Scalar(0.71007031138673404) * _tmp0 *
                         (-_tmp1 * _tmp32 -
                          std::sqrt(Scalar(std::pow(Scalar(-_tmp21 + p_b(1, 0)), Scalar(2)) +
                                           std::pow(Scalar(-_tmp31 + p_b(0, 0)), Scalar(2)))))));
  _res(1, 0) =
      -_tmp158 * (-_tmp129 * p_d(2, 0) + Scalar(1.0) * _tmp157 * std::sinh(_tmp155) -
                  (-_tmp129 * _tmp159 + _tmp161 * (-_tmp154 * _tmp160 - _tmp157 * _tmp158)) *
                      std::sinh(_tmp162)) -
      _tmp160 * (_tmp161 * p_d(2, 0) + std::cosh(_tmp155) - std::cosh(_tmp162));
  _res(2, 0) =
      -_tmp172 * (Scalar(1.0) * _tmp180 * std::sinh(_tmp170) - _tmp181 * p_a(2, 0) -
                  (_tmp171 * (-_tmp169 * _tmp176 - _tmp172 * _tmp180) - _tmp173 * _tmp181) *
                      std::sinh(_tmp174)) -
      _tmp176 * (_tmp171 * p_a(2, 0) + std::cosh(_tmp170) - std::cosh(_tmp174));
  _res(3, 0) = _tmp175 * (_tmp192 * p_c(2, 0) - std::cosh(_tmp193) + std::cosh(_tmp194)) -
               _tmp189 * (_tmp185 * p_c(2, 0) + Scalar(1.0) * _tmp191 * std::sinh(_tmp194) -
                          (_tmp185 * _tmp190 + _tmp192 * (_tmp175 * _tmp188 - _tmp189 * _tmp191)) *
                              std::sinh(_tmp193));

  return _res;
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym

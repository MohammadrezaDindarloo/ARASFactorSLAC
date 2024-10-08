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
 * Symbolic function: IK_residual_func_cost1_wrt_fh1_Nl22
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
Eigen::Matrix<Scalar, 4, 1> IkResidualFuncCost1WrtFh1Nl22(
    const Scalar fh1, const Scalar fv1, const sym::Rot3<Scalar>& DeltaRot,
    const Eigen::Matrix<Scalar, 3, 1>& position_vector, const Eigen::Matrix<Scalar, 4, 1>& encoder,
    const Eigen::Matrix<Scalar, 3, 1>& p_a, const Eigen::Matrix<Scalar, 3, 1>& p_b,
    const Eigen::Matrix<Scalar, 3, 1>& p_c, const Eigen::Matrix<Scalar, 3, 1>& p_d,
    const sym::Rot3<Scalar>& Rot_init, const Scalar epsilon) {
  // Total ops: 636

  // Unused inputs
  (void)encoder;
  (void)epsilon;

  // Input arrays
  const Eigen::Matrix<Scalar, 4, 1>& _DeltaRot = DeltaRot.Data();
  const Eigen::Matrix<Scalar, 4, 1>& _Rot_init = Rot_init.Data();

  // Intermediate terms (214)
  const Scalar _tmp0 = std::pow(fh1, Scalar(-2));
  const Scalar _tmp1 = Scalar(0.71007031138673404) * _tmp0;
  const Scalar _tmp2 = _DeltaRot[0] * _Rot_init[3] - _DeltaRot[1] * _Rot_init[2] +
                       _DeltaRot[2] * _Rot_init[1] + _DeltaRot[3] * _Rot_init[0];
  const Scalar _tmp3 = -2 * std::pow(_tmp2, Scalar(2));
  const Scalar _tmp4 = -_DeltaRot[0] * _Rot_init[1] + _DeltaRot[1] * _Rot_init[0] +
                       _DeltaRot[2] * _Rot_init[3] + _DeltaRot[3] * _Rot_init[2];
  const Scalar _tmp5 = 1 - 2 * std::pow(_tmp4, Scalar(2));
  const Scalar _tmp6 = Scalar(0.20999999999999999) * _tmp3 + Scalar(0.20999999999999999) * _tmp5;
  const Scalar _tmp7 = _DeltaRot[0] * _Rot_init[2] + _DeltaRot[1] * _Rot_init[3] -
                       _DeltaRot[2] * _Rot_init[0] + _DeltaRot[3] * _Rot_init[1];
  const Scalar _tmp8 = 2 * _tmp2 * _tmp7;
  const Scalar _tmp9 = -2 * _DeltaRot[0] * _Rot_init[0] - 2 * _DeltaRot[1] * _Rot_init[1] -
                       2 * _DeltaRot[2] * _Rot_init[2] + 2 * _DeltaRot[3] * _Rot_init[3];
  const Scalar _tmp10 = _tmp4 * _tmp9;
  const Scalar _tmp11 = Scalar(0.20999999999999999) * _tmp10 + Scalar(0.20999999999999999) * _tmp8;
  const Scalar _tmp12 = 2 * _tmp4;
  const Scalar _tmp13 = _tmp12 * _tmp7;
  const Scalar _tmp14 = _tmp2 * _tmp9;
  const Scalar _tmp15 = _tmp13 - _tmp14;
  const Scalar _tmp16 = -Scalar(0.010999999999999999) * _tmp15;
  const Scalar _tmp17 = -_tmp11 + _tmp16;
  const Scalar _tmp18 = _tmp17 + _tmp6;
  const Scalar _tmp19 = _tmp18 + position_vector(1, 0);
  const Scalar _tmp20 = -Scalar(0.20999999999999999) * _tmp10 + Scalar(0.20999999999999999) * _tmp8;
  const Scalar _tmp21 = -2 * std::pow(_tmp7, Scalar(2));
  const Scalar _tmp22 = Scalar(0.20999999999999999) * _tmp21 + Scalar(0.20999999999999999) * _tmp5;
  const Scalar _tmp23 = _tmp12 * _tmp2;
  const Scalar _tmp24 = _tmp7 * _tmp9;
  const Scalar _tmp25 = _tmp23 + _tmp24;
  const Scalar _tmp26 = -Scalar(0.010999999999999999) * _tmp25;
  const Scalar _tmp27 = -_tmp22 + _tmp26;
  const Scalar _tmp28 = _tmp20 + _tmp27;
  const Scalar _tmp29 = _tmp28 + position_vector(0, 0);
  const Scalar _tmp30 = Scalar(1.0) / (fh1);
  const Scalar _tmp31 = _tmp30 * fv1;
  const Scalar _tmp32 = std::asinh(_tmp31);
  const Scalar _tmp33 = Scalar(1.4083112389913199) * _tmp32;
  const Scalar _tmp34 =
      -_tmp33 * fh1 - std::sqrt(Scalar(std::pow(Scalar(-_tmp19 + p_d(1, 0)), Scalar(2)) +
                                       std::pow(Scalar(-_tmp29 + p_d(0, 0)), Scalar(2))));
  const Scalar _tmp35 = Scalar(0.71007031138673404) * _tmp30;
  const Scalar _tmp36 = _tmp34 * _tmp35;
  const Scalar _tmp37 =
      std::pow(Scalar(_tmp0 * std::pow(fv1, Scalar(2)) + 1), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp38 = Scalar(1.0) * _tmp32;
  const Scalar _tmp39 = Scalar(333.54000000000002) - fv1;
  const Scalar _tmp40 = Scalar(0.20999999999999999) * _tmp23 - Scalar(0.20999999999999999) * _tmp24;
  const Scalar _tmp41 = -Scalar(0.010999999999999999) * _tmp21 -
                        Scalar(0.010999999999999999) * _tmp3 + Scalar(-0.010999999999999999);
  const Scalar _tmp42 = Scalar(0.20999999999999999) * _tmp13 + Scalar(0.20999999999999999) * _tmp14;
  const Scalar _tmp43 = _tmp41 - _tmp42;
  const Scalar _tmp44 = _tmp40 + _tmp43;
  const Scalar _tmp45 = -_tmp20;
  const Scalar _tmp46 = _tmp22 + _tmp26;
  const Scalar _tmp47 = _tmp45 + _tmp46;
  const Scalar _tmp48 = _tmp47 + position_vector(0, 0);
  const Scalar _tmp49 = _tmp48 - p_b(0, 0);
  const Scalar _tmp50 = -_tmp6;
  const Scalar _tmp51 = _tmp11 + _tmp16;
  const Scalar _tmp52 = _tmp50 + _tmp51;
  const Scalar _tmp53 = _tmp52 + position_vector(1, 0);
  const Scalar _tmp54 = _tmp53 - p_b(1, 0);
  const Scalar _tmp55 = std::pow(Scalar(std::pow(_tmp49, Scalar(2)) + std::pow(_tmp54, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp56 = _tmp49 * _tmp55;
  const Scalar _tmp57 = _tmp54 * _tmp55;
  const Scalar _tmp58 = _tmp20 + _tmp46;
  const Scalar _tmp59 = _tmp58 + position_vector(0, 0);
  const Scalar _tmp60 = _tmp59 - p_c(0, 0);
  const Scalar _tmp61 = Scalar(1.0) / (_tmp60);
  const Scalar _tmp62 = _tmp51 + _tmp6;
  const Scalar _tmp63 = _tmp62 + position_vector(1, 0);
  const Scalar _tmp64 = _tmp63 - p_c(1, 0);
  const Scalar _tmp65 = _tmp61 * _tmp64;
  const Scalar _tmp66 = _tmp56 * _tmp65 - _tmp57;
  const Scalar _tmp67 = _tmp17 + _tmp50;
  const Scalar _tmp68 = _tmp67 + position_vector(1, 0);
  const Scalar _tmp69 = _tmp68 - p_a(1, 0);
  const Scalar _tmp70 = _tmp27 + _tmp45;
  const Scalar _tmp71 = _tmp70 + position_vector(0, 0);
  const Scalar _tmp72 = _tmp71 - p_a(0, 0);
  const Scalar _tmp73 = std::pow(Scalar(std::pow(_tmp69, Scalar(2)) + std::pow(_tmp72, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp74 = _tmp69 * _tmp73;
  const Scalar _tmp75 = _tmp72 * _tmp73;
  const Scalar _tmp76 = Scalar(1.0) / (_tmp65 * _tmp75 - _tmp74);
  const Scalar _tmp77 = -_tmp40;
  const Scalar _tmp78 = _tmp43 + _tmp77;
  const Scalar _tmp79 = _tmp41 + _tmp42;
  const Scalar _tmp80 = _tmp40 + _tmp79;
  const Scalar _tmp81 = _tmp76 * (-_tmp75 * _tmp78 + _tmp75 * _tmp80);
  const Scalar _tmp82 = _tmp56 * _tmp80;
  const Scalar _tmp83 = _tmp65 * _tmp80;
  const Scalar _tmp84 = _tmp76 * (_tmp74 * _tmp78 - _tmp75 * _tmp83);
  const Scalar _tmp85 = _tmp44 * _tmp57 - _tmp65 * _tmp82 - _tmp66 * _tmp84;
  const Scalar _tmp86 = Scalar(1.0) * _tmp62;
  const Scalar _tmp87 = -_tmp86;
  const Scalar _tmp88 = Scalar(1.0) / (_tmp67 + _tmp87);
  const Scalar _tmp89 = Scalar(1.0) * _tmp58;
  const Scalar _tmp90 = -_tmp70 + _tmp89;
  const Scalar _tmp91 = _tmp88 * _tmp90;
  const Scalar _tmp92 = -_tmp44 * _tmp56 - _tmp66 * _tmp81 + _tmp82 - _tmp85 * _tmp91;
  const Scalar _tmp93 = Scalar(1.0) / (_tmp92);
  const Scalar _tmp94 = _tmp86 * _tmp91 + _tmp89;
  const Scalar _tmp95 = 0;
  const Scalar _tmp96 = _tmp66 * _tmp76;
  const Scalar _tmp97 = _tmp75 * _tmp96;
  const Scalar _tmp98 =
      std::sqrt(Scalar(std::pow(_tmp60, Scalar(2)) + std::pow(_tmp64, Scalar(2))));
  const Scalar _tmp99 = _tmp61 * _tmp98;
  const Scalar _tmp100 = Scalar(1.0) / (_tmp98);
  const Scalar _tmp101 = _tmp99 * (_tmp100 * _tmp58 * _tmp64 - _tmp100 * _tmp60 * _tmp62);
  const Scalar _tmp102 = _tmp76 * (_tmp101 * _tmp75 + _tmp67 * _tmp75 - _tmp70 * _tmp74);
  const Scalar _tmp103 = Scalar(1.0) * _tmp88;
  const Scalar _tmp104 = _tmp103 * _tmp84 * _tmp90 - Scalar(1.0) * _tmp81;
  const Scalar _tmp105 = _tmp101 * _tmp56 - _tmp102 * _tmp66 - _tmp47 * _tmp57 + _tmp52 * _tmp56;
  const Scalar _tmp106 = _tmp105 * _tmp93;
  const Scalar _tmp107 = Scalar(1.0) / (_tmp105);
  const Scalar _tmp108 = _tmp107 * _tmp92;
  const Scalar _tmp109 = _tmp108 * (-Scalar(1.0) * _tmp102 - _tmp104 * _tmp106);
  const Scalar _tmp110 = _tmp93 * (_tmp104 + _tmp109);
  const Scalar _tmp111 = _tmp76 * (-_tmp110 * _tmp66 + Scalar(1.0));
  const Scalar _tmp112 = _tmp29 - p_d(0, 0);
  const Scalar _tmp113 = _tmp19 - p_d(1, 0);
  const Scalar _tmp114 =
      std::pow(Scalar(std::pow(_tmp112, Scalar(2)) + std::pow(_tmp113, Scalar(2))),
               Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp115 = _tmp113 * _tmp114;
  const Scalar _tmp116 = _tmp115 * _tmp99 * (_tmp110 * _tmp56 + _tmp111 * _tmp75);
  const Scalar _tmp117 = _tmp112 * _tmp114;
  const Scalar _tmp118 = _tmp115 * _tmp28 - _tmp117 * _tmp18;
  const Scalar _tmp119 = Scalar(1.0) * _tmp107;
  const Scalar _tmp120 = _tmp118 * _tmp99 * (_tmp119 * _tmp56 - _tmp119 * _tmp97);
  const Scalar _tmp121 = _tmp65 * _tmp84 + _tmp83;
  const Scalar _tmp122 = -_tmp121 * _tmp91 + _tmp65 * _tmp81 - _tmp80;
  const Scalar _tmp123 = _tmp108 * (-_tmp101 + _tmp102 * _tmp65 - _tmp106 * _tmp122);
  const Scalar _tmp124 = _tmp93 * (_tmp122 + _tmp123);
  const Scalar _tmp125 = _tmp76 * (-_tmp124 * _tmp66 - _tmp65);
  const Scalar _tmp126 = _tmp117 * _tmp99 * (_tmp124 * _tmp56 + _tmp125 * _tmp75 + Scalar(1.0));
  const Scalar _tmp127 = -_tmp116 * fh1 - _tmp120 * fh1 - _tmp126 * fh1 -
                         _tmp39 * _tmp99 * (_tmp56 * _tmp95 - _tmp95 * _tmp97);
  const Scalar _tmp128 = Scalar(1.0) / (_tmp127);
  const Scalar _tmp129 = _tmp52 + _tmp87;
  const Scalar _tmp130 = _tmp129 * _tmp91;
  const Scalar _tmp131 = Scalar(1.0) / (-_tmp130 - _tmp47 + _tmp89);
  const Scalar _tmp132 = Scalar(1.0) * _tmp131;
  const Scalar _tmp133 = _tmp129 * _tmp131;
  const Scalar _tmp134 = _tmp109 * _tmp133 - _tmp110 * _tmp85 - Scalar(1.0) * _tmp84;
  const Scalar _tmp135 = Scalar(1.0) * _tmp115 * (-_tmp103 * _tmp134 + _tmp109 * _tmp132);
  const Scalar _tmp136 = _tmp77 + _tmp79;
  const Scalar _tmp137 = _tmp136 * fh1;
  const Scalar _tmp138 = -_tmp115 * _tmp137 - Scalar(40.024799999999999) * _tmp15 - _tmp18 * fv1;
  const Scalar _tmp139 = _tmp130 * _tmp132 + Scalar(1.0);
  const Scalar _tmp140 = _tmp132 * _tmp91;
  const Scalar _tmp141 = -Scalar(1.0) * _tmp103 * _tmp139 + Scalar(1.0) * _tmp140;
  const Scalar _tmp142 = _tmp121 + _tmp123 * _tmp133 - _tmp124 * _tmp85;
  const Scalar _tmp143 = Scalar(1.0) * _tmp117 * (-_tmp103 * _tmp142 + _tmp123 * _tmp132);
  const Scalar _tmp144 = _tmp108 * _tmp132;
  const Scalar _tmp145 = -_tmp119 * _tmp85 + _tmp129 * _tmp144;
  const Scalar _tmp146 = Scalar(1.0) * _tmp118;
  const Scalar _tmp147 = _tmp146 * (-_tmp103 * _tmp145 + _tmp144);
  const Scalar _tmp148 = _tmp131 * _tmp94;
  const Scalar _tmp149 = -_tmp129 * _tmp148 - _tmp85 * _tmp95 + _tmp87;
  const Scalar _tmp150 = _tmp117 * _tmp137 + Scalar(40.024799999999999) * _tmp25 + _tmp28 * fv1;
  const Scalar _tmp151 = _tmp129 * _tmp88;
  const Scalar _tmp152 = Scalar(1.0) * _tmp132 * _tmp151 - Scalar(1.0) * _tmp132;
  const Scalar _tmp153 =
      _tmp135 * fh1 + _tmp138 * _tmp141 + _tmp143 * fh1 + _tmp147 * fh1 + _tmp150 * _tmp152 +
      Scalar(1.0) * _tmp39 * (-_tmp103 * _tmp149 - _tmp132 * _tmp94 + Scalar(1.0));
  const Scalar _tmp154 = std::asinh(_tmp128 * _tmp153);
  const Scalar _tmp155 = Scalar(1.4083112389913199) * _tmp127;
  const Scalar _tmp156 =
      -_tmp154 * _tmp155 - std::sqrt(Scalar(std::pow(Scalar(-_tmp59 + p_c(0, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp63 + p_c(1, 0)), Scalar(2))));
  const Scalar _tmp157 = Scalar(0.71007031138673404) * _tmp128;
  const Scalar _tmp158 = _tmp156 * _tmp157;
  const Scalar _tmp159 = Scalar(1.0) * _tmp154;
  const Scalar _tmp160 = -_tmp116 - _tmp120 - _tmp126;
  const Scalar _tmp161 = Scalar(1.4083112389913199) * _tmp160;
  const Scalar _tmp162 = std::pow(_tmp127, Scalar(-2));
  const Scalar _tmp163 = _tmp160 * _tmp162;
  const Scalar _tmp164 = _tmp115 * _tmp136;
  const Scalar _tmp165 = _tmp117 * _tmp136;
  const Scalar _tmp166 =
      (_tmp128 * (_tmp135 - _tmp141 * _tmp164 + _tmp143 + _tmp147 + _tmp152 * _tmp165) -
       _tmp153 * _tmp163) /
      std::sqrt(Scalar(std::pow(_tmp153, Scalar(2)) * _tmp162 + 1));
  const Scalar _tmp167 = Scalar(0.71007031138673404) * _tmp163;
  const Scalar _tmp168 = _tmp115 * _tmp134 * _tmp88;
  const Scalar _tmp169 = _tmp117 * _tmp88;
  const Scalar _tmp170 = _tmp142 * _tmp169;
  const Scalar _tmp171 = _tmp139 * _tmp88;
  const Scalar _tmp172 = _tmp118 * _tmp145 * _tmp88;
  const Scalar _tmp173 = _tmp132 * _tmp150;
  const Scalar _tmp174 = _tmp138 * _tmp171 + _tmp149 * _tmp39 * _tmp88 - _tmp151 * _tmp173 +
                         _tmp168 * fh1 + _tmp170 * fh1 + _tmp172 * fh1;
  const Scalar _tmp175 = _tmp117 * _tmp125;
  const Scalar _tmp176 = _tmp107 * _tmp146;
  const Scalar _tmp177 = _tmp176 * fh1;
  const Scalar _tmp178 = _tmp111 * _tmp115;
  const Scalar _tmp179 = _tmp39 * _tmp95;
  const Scalar _tmp180 = _tmp175 * fh1 - _tmp177 * _tmp96 + _tmp178 * fh1 - _tmp179 * _tmp96;
  const Scalar _tmp181 = Scalar(1.0) / (_tmp180);
  const Scalar _tmp182 = std::asinh(_tmp174 * _tmp181);
  const Scalar _tmp183 = Scalar(1.0) * _tmp182;
  const Scalar _tmp184 = std::pow(_tmp180, Scalar(-2));
  const Scalar _tmp185 = _tmp175 - _tmp176 * _tmp96 + _tmp178;
  const Scalar _tmp186 = _tmp184 * _tmp185;
  const Scalar _tmp187 =
      (-_tmp174 * _tmp186 + _tmp181 * (-_tmp129 * _tmp132 * _tmp136 * _tmp169 - _tmp164 * _tmp171 +
                                       _tmp168 + _tmp170 + _tmp172)) /
      std::sqrt(Scalar(std::pow(_tmp174, Scalar(2)) * _tmp184 + 1));
  const Scalar _tmp188 = Scalar(1.4083112389913199) * _tmp182;
  const Scalar _tmp189 =
      -_tmp180 * _tmp188 - std::sqrt(Scalar(std::pow(Scalar(-_tmp68 + p_a(1, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp71 + p_a(0, 0)), Scalar(2))));
  const Scalar _tmp190 = Scalar(0.71007031138673404) * _tmp181;
  const Scalar _tmp191 = _tmp189 * _tmp190;
  const Scalar _tmp192 = Scalar(1.4083112389913199) * _tmp180;
  const Scalar _tmp193 = Scalar(0.71007031138673404) * _tmp186;
  const Scalar _tmp194 = _tmp110 * _tmp115;
  const Scalar _tmp195 = _tmp117 * _tmp124;
  const Scalar _tmp196 = _tmp177 + _tmp179 + _tmp194 * fh1 + _tmp195 * fh1;
  const Scalar _tmp197 = std::pow(_tmp196, Scalar(-2));
  const Scalar _tmp198 = _tmp176 + _tmp194 + _tmp195;
  const Scalar _tmp199 = _tmp197 * _tmp198;
  const Scalar _tmp200 = Scalar(0.71007031138673404) * _tmp199;
  const Scalar _tmp201 = Scalar(1.0) / (_tmp196);
  const Scalar _tmp202 = _tmp118 * _tmp144;
  const Scalar _tmp203 = _tmp109 * _tmp115 * _tmp131;
  const Scalar _tmp204 = _tmp117 * _tmp123 * _tmp131;
  const Scalar _tmp205 = -_tmp138 * _tmp140 + _tmp148 * _tmp39 + _tmp173 - _tmp202 * fh1 -
                         _tmp203 * fh1 - _tmp204 * fh1;
  const Scalar _tmp206 = std::asinh(_tmp201 * _tmp205);
  const Scalar _tmp207 = Scalar(1.0) * _tmp206;
  const Scalar _tmp208 = (-_tmp199 * _tmp205 + _tmp201 * (_tmp132 * _tmp165 + _tmp140 * _tmp164 -
                                                          _tmp202 - _tmp203 - _tmp204)) /
                         std::sqrt(Scalar(_tmp197 * std::pow(_tmp205, Scalar(2)) + 1));
  const Scalar _tmp209 = Scalar(1.4083112389913199) * _tmp196;
  const Scalar _tmp210 =
      -_tmp206 * _tmp209 - std::sqrt(Scalar(std::pow(Scalar(-_tmp48 + p_b(0, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp53 + p_b(1, 0)), Scalar(2))));
  const Scalar _tmp211 = Scalar(1.4083112389913199) * _tmp198;
  const Scalar _tmp212 = Scalar(0.71007031138673404) * _tmp201;
  const Scalar _tmp213 = _tmp210 * _tmp212;

  // Output terms (1)
  Eigen::Matrix<Scalar, 4, 1> _res;

  _res(0, 0) =
      -Scalar(1.4083112389913199) * _tmp35 * p_d(2, 0) -
      Scalar(1.4083112389913199) * fh1 *
          (-Scalar(1.0) * _tmp0 * _tmp37 * fv1 * std::sinh(_tmp38) - _tmp1 * p_d(2, 0) -
           (-_tmp1 * _tmp34 + _tmp35 * (Scalar(1.4083112389913199) * _tmp31 * _tmp37 - _tmp33)) *
               std::sinh(_tmp36)) +
      Scalar(1.4083112389913199) * std::cosh(_tmp36) -
      Scalar(1.4083112389913199) * std::cosh(_tmp38);
  _res(1, 0) =
      -_tmp155 * (Scalar(1.0) * _tmp166 * std::sinh(_tmp159) - _tmp167 * p_c(2, 0) -
                  (-_tmp156 * _tmp167 + _tmp157 * (-_tmp154 * _tmp161 - _tmp155 * _tmp166)) *
                      std::sinh(_tmp158)) -
      _tmp161 * (_tmp157 * p_c(2, 0) - std::cosh(_tmp158) + std::cosh(_tmp159));
  _res(2, 0) =
      -Scalar(1.4083112389913199) * _tmp185 *
          (_tmp190 * p_a(2, 0) + std::cosh(_tmp183) - std::cosh(_tmp191)) -
      _tmp192 * (Scalar(1.0) * _tmp187 * std::sinh(_tmp183) - _tmp193 * p_a(2, 0) -
                 (-_tmp189 * _tmp193 + _tmp190 * (-_tmp185 * _tmp188 - _tmp187 * _tmp192)) *
                     std::sinh(_tmp191));
  _res(3, 0) =
      -_tmp209 * (-_tmp200 * p_b(2, 0) + Scalar(1.0) * _tmp208 * std::sinh(_tmp207) -
                  (-_tmp200 * _tmp210 + _tmp212 * (-_tmp206 * _tmp211 - _tmp208 * _tmp209)) *
                      std::sinh(_tmp213)) -
      _tmp211 * (_tmp212 * p_b(2, 0) + std::cosh(_tmp207) - std::cosh(_tmp213));

  return _res;
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym

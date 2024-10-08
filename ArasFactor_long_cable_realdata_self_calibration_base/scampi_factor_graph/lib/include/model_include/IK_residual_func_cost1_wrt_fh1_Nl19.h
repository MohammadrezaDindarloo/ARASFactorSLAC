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
 * Symbolic function: IK_residual_func_cost1_wrt_fh1_Nl19
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
Eigen::Matrix<Scalar, 4, 1> IkResidualFuncCost1WrtFh1Nl19(
    const Scalar fh1, const Scalar fv1, const sym::Rot3<Scalar>& DeltaRot,
    const Eigen::Matrix<Scalar, 3, 1>& position_vector, const Eigen::Matrix<Scalar, 4, 1>& encoder,
    const Eigen::Matrix<Scalar, 3, 1>& p_a, const Eigen::Matrix<Scalar, 3, 1>& p_b,
    const Eigen::Matrix<Scalar, 3, 1>& p_c, const Eigen::Matrix<Scalar, 3, 1>& p_d,
    const sym::Rot3<Scalar>& Rot_init, const Scalar epsilon) {
  // Total ops: 638

  // Unused inputs
  (void)encoder;
  (void)epsilon;

  // Input arrays
  const Eigen::Matrix<Scalar, 4, 1>& _DeltaRot = DeltaRot.Data();
  const Eigen::Matrix<Scalar, 4, 1>& _Rot_init = Rot_init.Data();

  // Intermediate terms (221)
  const Scalar _tmp0 = std::pow(fh1, Scalar(-2));
  const Scalar _tmp1 = Scalar(0.71007031138673404) * _tmp0;
  const Scalar _tmp2 = _DeltaRot[0] * _Rot_init[3] - _DeltaRot[1] * _Rot_init[2] +
                       _DeltaRot[2] * _Rot_init[1] + _DeltaRot[3] * _Rot_init[0];
  const Scalar _tmp3 = -2 * std::pow(_tmp2, Scalar(2));
  const Scalar _tmp4 = -_DeltaRot[0] * _Rot_init[1] + _DeltaRot[1] * _Rot_init[0] +
                       _DeltaRot[2] * _Rot_init[3] + _DeltaRot[3] * _Rot_init[2];
  const Scalar _tmp5 = -2 * std::pow(_tmp4, Scalar(2));
  const Scalar _tmp6 = Scalar(0.20999999999999999) * _tmp3 + Scalar(0.20999999999999999) * _tmp5 +
                       Scalar(0.20999999999999999);
  const Scalar _tmp7 = _DeltaRot[0] * _Rot_init[2] + _DeltaRot[1] * _Rot_init[3] -
                       _DeltaRot[2] * _Rot_init[0] + _DeltaRot[3] * _Rot_init[1];
  const Scalar _tmp8 = 2 * _tmp2;
  const Scalar _tmp9 = _tmp7 * _tmp8;
  const Scalar _tmp10 = -2 * _DeltaRot[0] * _Rot_init[0] - 2 * _DeltaRot[1] * _Rot_init[1] -
                        2 * _DeltaRot[2] * _Rot_init[2] + 2 * _DeltaRot[3] * _Rot_init[3];
  const Scalar _tmp11 = _tmp10 * _tmp4;
  const Scalar _tmp12 = Scalar(0.20999999999999999) * _tmp11 + Scalar(0.20999999999999999) * _tmp9;
  const Scalar _tmp13 = 2 * _tmp4 * _tmp7;
  const Scalar _tmp14 = _tmp10 * _tmp2;
  const Scalar _tmp15 = _tmp13 - _tmp14;
  const Scalar _tmp16 = -Scalar(0.010999999999999999) * _tmp15;
  const Scalar _tmp17 = -_tmp12 + _tmp16;
  const Scalar _tmp18 = _tmp17 + _tmp6;
  const Scalar _tmp19 = _tmp18 + position_vector(1, 0);
  const Scalar _tmp20 = 1 - 2 * std::pow(_tmp7, Scalar(2));
  const Scalar _tmp21 = Scalar(0.20999999999999999) * _tmp20 + Scalar(0.20999999999999999) * _tmp5;
  const Scalar _tmp22 = -_tmp21;
  const Scalar _tmp23 = _tmp4 * _tmp8;
  const Scalar _tmp24 = _tmp10 * _tmp7;
  const Scalar _tmp25 = _tmp23 + _tmp24;
  const Scalar _tmp26 = -Scalar(0.010999999999999999) * _tmp25;
  const Scalar _tmp27 = -Scalar(0.20999999999999999) * _tmp11 + Scalar(0.20999999999999999) * _tmp9;
  const Scalar _tmp28 = _tmp26 + _tmp27;
  const Scalar _tmp29 = _tmp22 + _tmp28;
  const Scalar _tmp30 = _tmp29 + position_vector(0, 0);
  const Scalar _tmp31 = Scalar(1.0) / (fh1);
  const Scalar _tmp32 = _tmp31 * fv1;
  const Scalar _tmp33 = std::asinh(_tmp32);
  const Scalar _tmp34 = Scalar(1.4083112389913199) * _tmp33;
  const Scalar _tmp35 =
      -_tmp34 * fh1 - std::sqrt(Scalar(std::pow(Scalar(-_tmp19 + p_d(1, 0)), Scalar(2)) +
                                       std::pow(Scalar(-_tmp30 + p_d(0, 0)), Scalar(2))));
  const Scalar _tmp36 = Scalar(0.71007031138673404) * _tmp31;
  const Scalar _tmp37 = _tmp35 * _tmp36;
  const Scalar _tmp38 =
      std::pow(Scalar(_tmp0 * std::pow(fv1, Scalar(2)) + 1), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp39 = Scalar(1.0) * _tmp33;
  const Scalar _tmp40 = _tmp26 - _tmp27;
  const Scalar _tmp41 = _tmp21 + _tmp40;
  const Scalar _tmp42 = _tmp41 + position_vector(0, 0);
  const Scalar _tmp43 = _tmp42 - p_b(0, 0);
  const Scalar _tmp44 = -_tmp6;
  const Scalar _tmp45 = _tmp12 + _tmp16;
  const Scalar _tmp46 = _tmp44 + _tmp45;
  const Scalar _tmp47 = _tmp46 + position_vector(1, 0);
  const Scalar _tmp48 = _tmp47 - p_b(1, 0);
  const Scalar _tmp49 = std::pow(Scalar(std::pow(_tmp43, Scalar(2)) + std::pow(_tmp48, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp50 = _tmp43 * _tmp49;
  const Scalar _tmp51 = Scalar(0.20999999999999999) * _tmp13 + Scalar(0.20999999999999999) * _tmp14;
  const Scalar _tmp52 = -_tmp51;
  const Scalar _tmp53 =
      -Scalar(0.010999999999999999) * _tmp20 - Scalar(0.010999999999999999) * _tmp3;
  const Scalar _tmp54 = Scalar(0.20999999999999999) * _tmp23 - Scalar(0.20999999999999999) * _tmp24;
  const Scalar _tmp55 = _tmp53 + _tmp54;
  const Scalar _tmp56 = _tmp52 + _tmp55;
  const Scalar _tmp57 = _tmp48 * _tmp49;
  const Scalar _tmp58 = _tmp51 + _tmp55;
  const Scalar _tmp59 = _tmp45 + _tmp6;
  const Scalar _tmp60 = _tmp59 + position_vector(1, 0);
  const Scalar _tmp61 = _tmp60 - p_c(1, 0);
  const Scalar _tmp62 = _tmp21 + _tmp28;
  const Scalar _tmp63 = _tmp62 + position_vector(0, 0);
  const Scalar _tmp64 = _tmp63 - p_c(0, 0);
  const Scalar _tmp65 = std::pow(Scalar(std::pow(_tmp61, Scalar(2)) + std::pow(_tmp64, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp66 = _tmp61 * _tmp65;
  const Scalar _tmp67 = _tmp53 - _tmp54;
  const Scalar _tmp68 = _tmp52 + _tmp67;
  const Scalar _tmp69 = _tmp22 + _tmp40;
  const Scalar _tmp70 = _tmp69 + position_vector(0, 0);
  const Scalar _tmp71 = _tmp70 - p_a(0, 0);
  const Scalar _tmp72 = Scalar(1.0) / (_tmp71);
  const Scalar _tmp73 = _tmp17 + _tmp44;
  const Scalar _tmp74 = _tmp73 + position_vector(1, 0);
  const Scalar _tmp75 = _tmp74 - p_a(1, 0);
  const Scalar _tmp76 = _tmp72 * _tmp75;
  const Scalar _tmp77 = _tmp68 * _tmp76;
  const Scalar _tmp78 = _tmp64 * _tmp65;
  const Scalar _tmp79 = _tmp58 * _tmp66 - _tmp77 * _tmp78;
  const Scalar _tmp80 = _tmp50 * _tmp76 - _tmp57;
  const Scalar _tmp81 = Scalar(1.0) / (-_tmp66 + _tmp76 * _tmp78);
  const Scalar _tmp82 = _tmp80 * _tmp81;
  const Scalar _tmp83 = _tmp50 * _tmp68;
  const Scalar _tmp84 = _tmp56 * _tmp57 - _tmp76 * _tmp83 - _tmp79 * _tmp82;
  const Scalar _tmp85 = Scalar(1.0) * _tmp73;
  const Scalar _tmp86 = -_tmp85;
  const Scalar _tmp87 = Scalar(1.0) / (_tmp59 + _tmp86);
  const Scalar _tmp88 = Scalar(1.0) * _tmp69;
  const Scalar _tmp89 = -_tmp62 + _tmp88;
  const Scalar _tmp90 = _tmp87 * _tmp89;
  const Scalar _tmp91 = -_tmp58 * _tmp78 + _tmp68 * _tmp78;
  const Scalar _tmp92 = -_tmp50 * _tmp56 - _tmp82 * _tmp91 + _tmp83 - _tmp84 * _tmp90;
  const Scalar _tmp93 = Scalar(1.0) / (_tmp92);
  const Scalar _tmp94 =
      std::sqrt(Scalar(std::pow(_tmp71, Scalar(2)) + std::pow(_tmp75, Scalar(2))));
  const Scalar _tmp95 = Scalar(1.0) / (_tmp94);
  const Scalar _tmp96 = _tmp72 * _tmp94;
  const Scalar _tmp97 = _tmp96 * (_tmp69 * _tmp75 * _tmp95 - _tmp71 * _tmp73 * _tmp95);
  const Scalar _tmp98 = _tmp59 * _tmp78 - _tmp62 * _tmp66 + _tmp78 * _tmp97;
  const Scalar _tmp99 = _tmp76 * _tmp81;
  const Scalar _tmp100 = _tmp77 + _tmp79 * _tmp99;
  const Scalar _tmp101 = -_tmp100 * _tmp90 - _tmp68 + _tmp91 * _tmp99;
  const Scalar _tmp102 = -_tmp41 * _tmp57 + _tmp46 * _tmp50 + _tmp50 * _tmp97 - _tmp82 * _tmp98;
  const Scalar _tmp103 = _tmp102 * _tmp93;
  const Scalar _tmp104 = Scalar(1.0) / (_tmp102);
  const Scalar _tmp105 = _tmp104 * _tmp92;
  const Scalar _tmp106 = _tmp105 * (-_tmp101 * _tmp103 - _tmp97 + _tmp98 * _tmp99);
  const Scalar _tmp107 = _tmp101 + _tmp106;
  const Scalar _tmp108 = _tmp107 * _tmp93;
  const Scalar _tmp109 = _tmp80 * _tmp93;
  const Scalar _tmp110 = _tmp81 * (-_tmp107 * _tmp109 - _tmp76);
  const Scalar _tmp111 = _tmp30 - p_d(0, 0);
  const Scalar _tmp112 = _tmp19 - p_d(1, 0);
  const Scalar _tmp113 =
      std::pow(Scalar(std::pow(_tmp111, Scalar(2)) + std::pow(_tmp112, Scalar(2))),
               Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp114 = _tmp111 * _tmp113;
  const Scalar _tmp115 = _tmp114 * _tmp96 * (_tmp108 * _tmp50 + _tmp110 * _tmp78 + Scalar(1.0));
  const Scalar _tmp116 = Scalar(333.54000000000002) - fv1;
  const Scalar _tmp117 = _tmp85 * _tmp90 + _tmp88;
  const Scalar _tmp118 = 0;
  const Scalar _tmp119 = _tmp118 * _tmp93;
  const Scalar _tmp120 = _tmp78 * _tmp82;
  const Scalar _tmp121 = _tmp112 * _tmp113;
  const Scalar _tmp122 = -_tmp114 * _tmp18 + _tmp121 * _tmp29;
  const Scalar _tmp123 = Scalar(1.0) * _tmp104;
  const Scalar _tmp124 = _tmp122 * _tmp96 * (-_tmp120 * _tmp123 + _tmp123 * _tmp50);
  const Scalar _tmp125 = Scalar(1.0) * _tmp81;
  const Scalar _tmp126 = Scalar(1.0) * _tmp87;
  const Scalar _tmp127 = -_tmp125 * _tmp91 + _tmp126 * _tmp79 * _tmp81 * _tmp89;
  const Scalar _tmp128 = _tmp105 * (-_tmp103 * _tmp127 - _tmp125 * _tmp98);
  const Scalar _tmp129 = _tmp127 + _tmp128;
  const Scalar _tmp130 = _tmp81 * (-_tmp109 * _tmp129 + Scalar(1.0));
  const Scalar _tmp131 = _tmp129 * _tmp93;
  const Scalar _tmp132 = _tmp121 * _tmp96 * (_tmp130 * _tmp78 + _tmp131 * _tmp50);
  const Scalar _tmp133 = -_tmp115 * fh1 -
                         _tmp116 * _tmp96 * (-_tmp119 * _tmp120 + _tmp119 * _tmp50) -
                         _tmp124 * fh1 - _tmp132 * fh1;
  const Scalar _tmp134 = Scalar(1.0) / (_tmp133);
  const Scalar _tmp135 = _tmp51 + _tmp67;
  const Scalar _tmp136 = _tmp135 * fh1;
  const Scalar _tmp137 = -_tmp121 * _tmp136 - Scalar(40.024799999999999) * _tmp15 - _tmp18 * fv1;
  const Scalar _tmp138 = _tmp46 + _tmp86;
  const Scalar _tmp139 = _tmp138 * _tmp90;
  const Scalar _tmp140 = Scalar(1.0) / (-_tmp139 - _tmp41 + _tmp88);
  const Scalar _tmp141 = Scalar(1.0) * _tmp140;
  const Scalar _tmp142 = _tmp139 * _tmp141 + Scalar(1.0);
  const Scalar _tmp143 = _tmp141 * _tmp90;
  const Scalar _tmp144 = -Scalar(1.0) * _tmp126 * _tmp142 + Scalar(1.0) * _tmp143;
  const Scalar _tmp145 = _tmp105 * _tmp141;
  const Scalar _tmp146 = -_tmp123 * _tmp84 + _tmp138 * _tmp145;
  const Scalar _tmp147 = Scalar(1.0) * _tmp122;
  const Scalar _tmp148 = _tmp147 * (-_tmp126 * _tmp146 + _tmp145);
  const Scalar _tmp149 = _tmp138 * _tmp140;
  const Scalar _tmp150 = _tmp84 * _tmp93;
  const Scalar _tmp151 = -_tmp125 * _tmp79 + _tmp128 * _tmp149 - _tmp129 * _tmp150;
  const Scalar _tmp152 = Scalar(1.0) * _tmp121 * (-_tmp126 * _tmp151 + _tmp128 * _tmp141);
  const Scalar _tmp153 = _tmp100 + _tmp106 * _tmp149 - _tmp107 * _tmp150;
  const Scalar _tmp154 = Scalar(1.0) * _tmp114 * (_tmp106 * _tmp141 - _tmp126 * _tmp153);
  const Scalar _tmp155 = _tmp114 * _tmp136 + Scalar(40.024799999999999) * _tmp25 + _tmp29 * fv1;
  const Scalar _tmp156 = _tmp138 * _tmp87;
  const Scalar _tmp157 = Scalar(1.0) * _tmp141 * _tmp156 - Scalar(1.0) * _tmp141;
  const Scalar _tmp158 = _tmp117 * _tmp140;
  const Scalar _tmp159 = -_tmp118 * _tmp150 - _tmp138 * _tmp158 + _tmp86;
  const Scalar _tmp160 =
      Scalar(1.0) * _tmp116 * (-_tmp117 * _tmp141 - _tmp126 * _tmp159 + Scalar(1.0)) +
      _tmp137 * _tmp144 + _tmp148 * fh1 + _tmp152 * fh1 + _tmp154 * fh1 + _tmp155 * _tmp157;
  const Scalar _tmp161 = std::asinh(_tmp134 * _tmp160);
  const Scalar _tmp162 = Scalar(1.0) * _tmp161;
  const Scalar _tmp163 = _tmp121 * _tmp135;
  const Scalar _tmp164 = _tmp114 * _tmp135;
  const Scalar _tmp165 = std::pow(_tmp133, Scalar(-2));
  const Scalar _tmp166 = -_tmp115 - _tmp124 - _tmp132;
  const Scalar _tmp167 = _tmp165 * _tmp166;
  const Scalar _tmp168 =
      (_tmp134 * (-_tmp144 * _tmp163 + _tmp148 + _tmp152 + _tmp154 + _tmp157 * _tmp164) -
       _tmp160 * _tmp167) /
      std::sqrt(Scalar(std::pow(_tmp160, Scalar(2)) * _tmp165 + 1));
  const Scalar _tmp169 = Scalar(0.71007031138673404) * _tmp167;
  const Scalar _tmp170 = Scalar(1.4083112389913199) * _tmp133;
  const Scalar _tmp171 =
      -_tmp161 * _tmp170 - std::sqrt(Scalar(std::pow(Scalar(-_tmp70 + p_a(0, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp74 + p_a(1, 0)), Scalar(2))));
  const Scalar _tmp172 = Scalar(0.71007031138673404) * _tmp134;
  const Scalar _tmp173 = _tmp171 * _tmp172;
  const Scalar _tmp174 = Scalar(1.4083112389913199) * _tmp166;
  const Scalar _tmp175 = _tmp121 * _tmp130;
  const Scalar _tmp176 = _tmp104 * _tmp147;
  const Scalar _tmp177 = _tmp176 * fh1;
  const Scalar _tmp178 = _tmp110 * _tmp114;
  const Scalar _tmp179 = _tmp116 * _tmp119;
  const Scalar _tmp180 = _tmp175 * fh1 - _tmp177 * _tmp82 + _tmp178 * fh1 - _tmp179 * _tmp82;
  const Scalar _tmp181 = Scalar(1.0) / (_tmp180);
  const Scalar _tmp182 = _tmp122 * _tmp146 * _tmp87;
  const Scalar _tmp183 = _tmp142 * _tmp87;
  const Scalar _tmp184 = _tmp121 * _tmp151 * _tmp87;
  const Scalar _tmp185 = _tmp141 * _tmp155;
  const Scalar _tmp186 = _tmp114 * _tmp87;
  const Scalar _tmp187 = _tmp153 * _tmp186;
  const Scalar _tmp188 = _tmp116 * _tmp159 * _tmp87 + _tmp137 * _tmp183 - _tmp156 * _tmp185 +
                         _tmp182 * fh1 + _tmp184 * fh1 + _tmp187 * fh1;
  const Scalar _tmp189 = std::asinh(_tmp181 * _tmp188);
  const Scalar _tmp190 = Scalar(1.0) * _tmp189;
  const Scalar _tmp191 = Scalar(0.71007031138673404) * _tmp181;
  const Scalar _tmp192 = Scalar(1.4083112389913199) * _tmp180;
  const Scalar _tmp193 =
      -_tmp189 * _tmp192 - std::sqrt(Scalar(std::pow(Scalar(-_tmp60 + p_c(1, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp63 + p_c(0, 0)), Scalar(2))));
  const Scalar _tmp194 = _tmp191 * _tmp193;
  const Scalar _tmp195 = _tmp175 - _tmp176 * _tmp82 + _tmp178;
  const Scalar _tmp196 = Scalar(1.4083112389913199) * _tmp195;
  const Scalar _tmp197 = std::pow(_tmp180, Scalar(-2));
  const Scalar _tmp198 = _tmp195 * _tmp197;
  const Scalar _tmp199 = Scalar(0.71007031138673404) * _tmp198;
  const Scalar _tmp200 = (_tmp181 * (-_tmp135 * _tmp138 * _tmp141 * _tmp186 - _tmp163 * _tmp183 +
                                     _tmp182 + _tmp184 + _tmp187) -
                          _tmp188 * _tmp198) /
                         std::sqrt(Scalar(std::pow(_tmp188, Scalar(2)) * _tmp197 + 1));
  const Scalar _tmp201 = _tmp121 * _tmp131;
  const Scalar _tmp202 = _tmp108 * _tmp114;
  const Scalar _tmp203 = _tmp177 + _tmp179 + _tmp201 * fh1 + _tmp202 * fh1;
  const Scalar _tmp204 = Scalar(1.0) / (_tmp203);
  const Scalar _tmp205 = _tmp106 * _tmp114 * _tmp140;
  const Scalar _tmp206 = _tmp122 * _tmp145;
  const Scalar _tmp207 = _tmp121 * _tmp128 * _tmp140;
  const Scalar _tmp208 = _tmp116 * _tmp158 - _tmp137 * _tmp143 + _tmp185 - _tmp205 * fh1 -
                         _tmp206 * fh1 - _tmp207 * fh1;
  const Scalar _tmp209 = std::asinh(_tmp204 * _tmp208);
  const Scalar _tmp210 = Scalar(1.4083112389913199) * _tmp203;
  const Scalar _tmp211 =
      -_tmp209 * _tmp210 - std::sqrt(Scalar(std::pow(Scalar(-_tmp42 + p_b(0, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp47 + p_b(1, 0)), Scalar(2))));
  const Scalar _tmp212 = std::pow(_tmp203, Scalar(-2));
  const Scalar _tmp213 = _tmp176 + _tmp201 + _tmp202;
  const Scalar _tmp214 = _tmp212 * _tmp213;
  const Scalar _tmp215 = Scalar(0.71007031138673404) * _tmp214;
  const Scalar _tmp216 = Scalar(1.4083112389913199) * _tmp213;
  const Scalar _tmp217 =
      (_tmp204 * (_tmp141 * _tmp164 + _tmp143 * _tmp163 - _tmp205 - _tmp206 - _tmp207) -
       _tmp208 * _tmp214) /
      std::sqrt(Scalar(std::pow(_tmp208, Scalar(2)) * _tmp212 + 1));
  const Scalar _tmp218 = Scalar(0.71007031138673404) * _tmp204;
  const Scalar _tmp219 = _tmp211 * _tmp218;
  const Scalar _tmp220 = Scalar(1.0) * _tmp209;

  // Output terms (1)
  Eigen::Matrix<Scalar, 4, 1> _res;

  _res(0, 0) =
      -Scalar(1.4083112389913199) * _tmp36 * p_d(2, 0) -
      Scalar(1.4083112389913199) * fh1 *
          (-Scalar(1.0) * _tmp0 * _tmp38 * fv1 * std::sinh(_tmp39) - _tmp1 * p_d(2, 0) -
           (-_tmp1 * _tmp35 + _tmp36 * (Scalar(1.4083112389913199) * _tmp32 * _tmp38 - _tmp34)) *
               std::sinh(_tmp37)) +
      Scalar(1.4083112389913199) * std::cosh(_tmp37) -
      Scalar(1.4083112389913199) * std::cosh(_tmp39);
  _res(1, 0) =
      -_tmp170 * (Scalar(1.0) * _tmp168 * std::sinh(_tmp162) - _tmp169 * p_a(2, 0) -
                  (-_tmp169 * _tmp171 + _tmp172 * (-_tmp161 * _tmp174 - _tmp168 * _tmp170)) *
                      std::sinh(_tmp173)) -
      _tmp174 * (_tmp172 * p_a(2, 0) + std::cosh(_tmp162) - std::cosh(_tmp173));
  _res(2, 0) =
      -_tmp192 * (-_tmp199 * p_c(2, 0) + Scalar(1.0) * _tmp200 * std::sinh(_tmp190) -
                  (_tmp191 * (-_tmp189 * _tmp196 - _tmp192 * _tmp200) - _tmp193 * _tmp199) *
                      std::sinh(_tmp194)) -
      _tmp196 * (_tmp191 * p_c(2, 0) + std::cosh(_tmp190) - std::cosh(_tmp194));
  _res(3, 0) =
      -_tmp210 * (-_tmp215 * p_b(2, 0) + Scalar(1.0) * _tmp217 * std::sinh(_tmp220) -
                  (-_tmp211 * _tmp215 + _tmp218 * (-_tmp209 * _tmp216 - _tmp210 * _tmp217)) *
                      std::sinh(_tmp219)) -
      _tmp216 * (_tmp218 * p_b(2, 0) - std::cosh(_tmp219) + std::cosh(_tmp220));

  return _res;
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym

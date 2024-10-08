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
 * Symbolic function: IK_residual_func_cost1_wrt_pa_Nl12
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
 *     res: Matrix43
 */
template <typename Scalar>
Eigen::Matrix<Scalar, 4, 3> IkResidualFuncCost1WrtPaNl12(
    const Scalar fh1, const Scalar fv1, const sym::Rot3<Scalar>& DeltaRot,
    const Eigen::Matrix<Scalar, 3, 1>& position_vector, const Eigen::Matrix<Scalar, 4, 1>& encoder,
    const Eigen::Matrix<Scalar, 3, 1>& p_a, const Eigen::Matrix<Scalar, 3, 1>& p_b,
    const Eigen::Matrix<Scalar, 3, 1>& p_c, const Eigen::Matrix<Scalar, 3, 1>& p_d,
    const sym::Rot3<Scalar>& Rot_init, const Scalar epsilon) {
  // Total ops: 1219

  // Unused inputs
  (void)encoder;
  (void)epsilon;

  // Input arrays
  const Eigen::Matrix<Scalar, 4, 1>& _DeltaRot = DeltaRot.Data();
  const Eigen::Matrix<Scalar, 4, 1>& _Rot_init = Rot_init.Data();

  // Intermediate terms (395)
  const Scalar _tmp0 = Scalar(333.54000000000002) - fv1;
  const Scalar _tmp1 = _DeltaRot[0] * _Rot_init[2] + _DeltaRot[1] * _Rot_init[3] -
                       _DeltaRot[2] * _Rot_init[0] + _DeltaRot[3] * _Rot_init[1];
  const Scalar _tmp2 = -_DeltaRot[0] * _Rot_init[1] + _DeltaRot[1] * _Rot_init[0] +
                       _DeltaRot[2] * _Rot_init[3] + _DeltaRot[3] * _Rot_init[2];
  const Scalar _tmp3 = 2 * _tmp2;
  const Scalar _tmp4 = _tmp1 * _tmp3;
  const Scalar _tmp5 = _DeltaRot[0] * _Rot_init[3] - _DeltaRot[1] * _Rot_init[2] +
                       _DeltaRot[2] * _Rot_init[1] + _DeltaRot[3] * _Rot_init[0];
  const Scalar _tmp6 = -2 * _DeltaRot[0] * _Rot_init[0] - 2 * _DeltaRot[1] * _Rot_init[1] -
                       2 * _DeltaRot[2] * _Rot_init[2] + 2 * _DeltaRot[3] * _Rot_init[3];
  const Scalar _tmp7 = _tmp5 * _tmp6;
  const Scalar _tmp8 = Scalar(0.20999999999999999) * _tmp4 + Scalar(0.20999999999999999) * _tmp7;
  const Scalar _tmp9 = -_tmp8;
  const Scalar _tmp10 = -2 * std::pow(_tmp1, Scalar(2));
  const Scalar _tmp11 = -2 * std::pow(_tmp5, Scalar(2));
  const Scalar _tmp12 = -Scalar(0.010999999999999999) * _tmp10 -
                        Scalar(0.010999999999999999) * _tmp11 + Scalar(-0.010999999999999999);
  const Scalar _tmp13 = _tmp3 * _tmp5;
  const Scalar _tmp14 = _tmp1 * _tmp6;
  const Scalar _tmp15 = Scalar(0.20999999999999999) * _tmp13 - Scalar(0.20999999999999999) * _tmp14;
  const Scalar _tmp16 = _tmp12 - _tmp15;
  const Scalar _tmp17 = _tmp16 + _tmp9;
  const Scalar _tmp18 = 2 * _tmp1 * _tmp5;
  const Scalar _tmp19 = _tmp2 * _tmp6;
  const Scalar _tmp20 = Scalar(0.20999999999999999) * _tmp18 - Scalar(0.20999999999999999) * _tmp19;
  const Scalar _tmp21 = _tmp13 + _tmp14;
  const Scalar _tmp22 = -Scalar(0.010999999999999999) * _tmp21;
  const Scalar _tmp23 = 1 - 2 * std::pow(_tmp2, Scalar(2));
  const Scalar _tmp24 = Scalar(0.20999999999999999) * _tmp10 + Scalar(0.20999999999999999) * _tmp23;
  const Scalar _tmp25 = _tmp22 - _tmp24;
  const Scalar _tmp26 = _tmp20 + _tmp25;
  const Scalar _tmp27 = _tmp26 + position_vector(0, 0);
  const Scalar _tmp28 = _tmp27 - p_d(0, 0);
  const Scalar _tmp29 = Scalar(0.20999999999999999) * _tmp11 + Scalar(0.20999999999999999) * _tmp23;
  const Scalar _tmp30 = Scalar(0.20999999999999999) * _tmp18 + Scalar(0.20999999999999999) * _tmp19;
  const Scalar _tmp31 = _tmp4 - _tmp7;
  const Scalar _tmp32 = -Scalar(0.010999999999999999) * _tmp31;
  const Scalar _tmp33 = -_tmp30 + _tmp32;
  const Scalar _tmp34 = _tmp29 + _tmp33;
  const Scalar _tmp35 = _tmp34 + position_vector(1, 0);
  const Scalar _tmp36 = _tmp35 - p_d(1, 0);
  const Scalar _tmp37 = std::pow(Scalar(std::pow(_tmp28, Scalar(2)) + std::pow(_tmp36, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp38 = _tmp28 * _tmp37;
  const Scalar _tmp39 = _tmp17 * _tmp38;
  const Scalar _tmp40 = _tmp12 + _tmp15;
  const Scalar _tmp41 = _tmp40 + _tmp9;
  const Scalar _tmp42 = -_tmp20;
  const Scalar _tmp43 = _tmp22 + _tmp24;
  const Scalar _tmp44 = _tmp42 + _tmp43;
  const Scalar _tmp45 = _tmp44 + position_vector(0, 0);
  const Scalar _tmp46 = _tmp45 - p_b(0, 0);
  const Scalar _tmp47 = std::pow(_tmp46, Scalar(2));
  const Scalar _tmp48 = -_tmp29;
  const Scalar _tmp49 = _tmp30 + _tmp32;
  const Scalar _tmp50 = _tmp48 + _tmp49;
  const Scalar _tmp51 = _tmp50 + position_vector(1, 0);
  const Scalar _tmp52 = _tmp51 - p_b(1, 0);
  const Scalar _tmp53 = _tmp47 + std::pow(_tmp52, Scalar(2));
  const Scalar _tmp54 = std::pow(_tmp53, Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp55 = _tmp46 * _tmp54;
  const Scalar _tmp56 = _tmp17 * _tmp55;
  const Scalar _tmp57 = -_tmp41 * _tmp55 + _tmp56;
  const Scalar _tmp58 = _tmp52 * _tmp54;
  const Scalar _tmp59 = _tmp25 + _tmp42;
  const Scalar _tmp60 = _tmp59 + position_vector(0, 0);
  const Scalar _tmp61 = _tmp60 - p_a(0, 0);
  const Scalar _tmp62 = Scalar(1.0) / (_tmp61);
  const Scalar _tmp63 = _tmp33 + _tmp48;
  const Scalar _tmp64 = _tmp63 + position_vector(1, 0);
  const Scalar _tmp65 = _tmp64 - p_a(1, 0);
  const Scalar _tmp66 = _tmp62 * _tmp65;
  const Scalar _tmp67 = _tmp55 * _tmp66;
  const Scalar _tmp68 = -_tmp58 + _tmp67;
  const Scalar _tmp69 = Scalar(1.0) / (_tmp68);
  const Scalar _tmp70 = _tmp36 * _tmp37;
  const Scalar _tmp71 = _tmp38 * _tmp66;
  const Scalar _tmp72 = -_tmp70 + _tmp71;
  const Scalar _tmp73 = _tmp69 * _tmp72;
  const Scalar _tmp74 = _tmp16 + _tmp8;
  const Scalar _tmp75 = _tmp17 * _tmp62;
  const Scalar _tmp76 = _tmp65 * _tmp75;
  const Scalar _tmp77 = _tmp41 * _tmp58 - _tmp55 * _tmp76;
  const Scalar _tmp78 = _tmp69 * _tmp77;
  const Scalar _tmp79 = -_tmp38 * _tmp76 + _tmp70 * _tmp74 - _tmp72 * _tmp78;
  const Scalar _tmp80 = Scalar(1.0) * _tmp63;
  const Scalar _tmp81 = -_tmp80;
  const Scalar _tmp82 = Scalar(1.0) / (_tmp50 + _tmp81);
  const Scalar _tmp83 = Scalar(1.0) * _tmp59;
  const Scalar _tmp84 = _tmp82 * (-_tmp44 + _tmp83);
  const Scalar _tmp85 = -_tmp38 * _tmp74 + _tmp39 - _tmp57 * _tmp73 - _tmp79 * _tmp84;
  const Scalar _tmp86 = Scalar(1.0) / (_tmp85);
  const Scalar _tmp87 = _tmp80 * _tmp84 + _tmp83;
  const Scalar _tmp88 = 0;
  const Scalar _tmp89 = _tmp86 * _tmp88;
  const Scalar _tmp90 = _tmp38 * _tmp89;
  const Scalar _tmp91 = _tmp55 * _tmp73;
  const Scalar _tmp92 = _tmp0 * (-_tmp89 * _tmp91 + _tmp90);
  const Scalar _tmp93 = std::pow(_tmp61, Scalar(2));
  const Scalar _tmp94 = Scalar(1.0) / (_tmp93);
  const Scalar _tmp95 = std::pow(_tmp65, Scalar(2));
  const Scalar _tmp96 = _tmp93 + _tmp95;
  const Scalar _tmp97 = std::sqrt(_tmp96);
  const Scalar _tmp98 = _tmp94 * _tmp97;
  const Scalar _tmp99 = Scalar(1.0) / (_tmp97);
  const Scalar _tmp100 = _tmp63 * _tmp99;
  const Scalar _tmp101 = _tmp59 * _tmp99;
  const Scalar _tmp102 = -_tmp100 * _tmp61 + _tmp101 * _tmp65;
  const Scalar _tmp103 = _tmp62 * _tmp97;
  const Scalar _tmp104 = _tmp102 * _tmp103;
  const Scalar _tmp105 = _tmp104 * _tmp55 - _tmp44 * _tmp58 + _tmp50 * _tmp55;
  const Scalar _tmp106 = Scalar(1.0) * _tmp69;
  const Scalar _tmp107 = _tmp105 * _tmp69;
  const Scalar _tmp108 = _tmp104 * _tmp38 - _tmp107 * _tmp72 - _tmp26 * _tmp70 + _tmp34 * _tmp38;
  const Scalar _tmp109 = _tmp106 * _tmp77;
  const Scalar _tmp110 = -_tmp106 * _tmp57 + _tmp109 * _tmp84;
  const Scalar _tmp111 = _tmp110 * _tmp86;
  const Scalar _tmp112 = -_tmp105 * _tmp106 - _tmp108 * _tmp111;
  const Scalar _tmp113 = Scalar(1.0) / (_tmp108);
  const Scalar _tmp114 = _tmp113 * _tmp85;
  const Scalar _tmp115 = _tmp112 * _tmp114;
  const Scalar _tmp116 = _tmp110 + _tmp115;
  const Scalar _tmp117 = _tmp116 * _tmp86;
  const Scalar _tmp118 = -_tmp117 * _tmp72 + Scalar(1.0);
  const Scalar _tmp119 = _tmp55 * _tmp69;
  const Scalar _tmp120 = _tmp117 * _tmp38;
  const Scalar _tmp121 = _tmp29 + _tmp49;
  const Scalar _tmp122 = _tmp121 - p_c(1, 0) + position_vector(1, 0);
  const Scalar _tmp123 = _tmp20 + _tmp43;
  const Scalar _tmp124 = _tmp123 - p_c(0, 0) + position_vector(0, 0);
  const Scalar _tmp125 =
      std::pow(Scalar(std::pow(_tmp122, Scalar(2)) + std::pow(_tmp124, Scalar(2))),
               Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp126 = _tmp122 * _tmp125;
  const Scalar _tmp127 = _tmp126 * fh1;
  const Scalar _tmp128 = _tmp127 * (_tmp118 * _tmp119 + _tmp120);
  const Scalar _tmp129 = _tmp128 * _tmp99;
  const Scalar _tmp130 = _tmp106 * _tmp55;
  const Scalar _tmp131 = std::pow(_tmp96, Scalar(Scalar(-3) / Scalar(2)));
  const Scalar _tmp132 = _tmp131 * _tmp63;
  const Scalar _tmp133 = _tmp131 * _tmp59;
  const Scalar _tmp134 = _tmp61 * _tmp65;
  const Scalar _tmp135 = _tmp103 * (_tmp100 - _tmp132 * _tmp93 + _tmp133 * _tmp134);
  const Scalar _tmp136 = _tmp102 * _tmp99;
  const Scalar _tmp137 = _tmp102 * _tmp98;
  const Scalar _tmp138 = _tmp135 * _tmp55 - _tmp136 * _tmp55 + _tmp137 * _tmp55;
  const Scalar _tmp139 = std::pow(_tmp68, Scalar(-2));
  const Scalar _tmp140 = _tmp65 * _tmp94;
  const Scalar _tmp141 = _tmp139 * _tmp140;
  const Scalar _tmp142 = _tmp105 * _tmp55;
  const Scalar _tmp143 = _tmp141 * _tmp142;
  const Scalar _tmp144 = _tmp107 * _tmp140;
  const Scalar _tmp145 = _tmp135 * _tmp38 - _tmp136 * _tmp38 + _tmp137 * _tmp38 - _tmp138 * _tmp73 +
                         _tmp143 * _tmp72 - _tmp144 * _tmp38;
  const Scalar _tmp146 = std::pow(_tmp108, Scalar(-2));
  const Scalar _tmp147 = _tmp145 * _tmp146;
  const Scalar _tmp148 = _tmp147 * _tmp72;
  const Scalar _tmp149 = Scalar(1.0) * _tmp113;
  const Scalar _tmp150 = _tmp47 / _tmp53;
  const Scalar _tmp151 = _tmp141 * _tmp150;
  const Scalar _tmp152 = _tmp151 * _tmp72;
  const Scalar _tmp153 = Scalar(1.0) * _tmp38;
  const Scalar _tmp154 = _tmp149 * _tmp38;
  const Scalar _tmp155 = _tmp119 * _tmp140;
  const Scalar _tmp156 = _tmp124 * _tmp125;
  const Scalar _tmp157 = fh1 * (-_tmp121 * _tmp156 + _tmp123 * _tmp126);
  const Scalar _tmp158 = _tmp103 * _tmp157;
  const Scalar _tmp159 = _tmp66 * _tmp78 + _tmp76;
  const Scalar _tmp160 = _tmp57 * _tmp69;
  const Scalar _tmp161 = -_tmp159 * _tmp84 + _tmp160 * _tmp66 - _tmp17;
  const Scalar _tmp162 = _tmp161 * _tmp86;
  const Scalar _tmp163 = -_tmp104 + _tmp107 * _tmp66 - _tmp108 * _tmp162;
  const Scalar _tmp164 = _tmp113 * _tmp163;
  const Scalar _tmp165 = _tmp164 * _tmp85;
  const Scalar _tmp166 = _tmp161 + _tmp165;
  const Scalar _tmp167 = _tmp72 * _tmp86;
  const Scalar _tmp168 = -_tmp166 * _tmp167 - _tmp66;
  const Scalar _tmp169 = _tmp95 / [&]() {
    const Scalar base = _tmp61;
    return base * base * base;
  }();
  const Scalar _tmp170 = _tmp139 * _tmp169;
  const Scalar _tmp171 = _tmp55 * _tmp77;
  const Scalar _tmp172 = _tmp141 * _tmp171;
  const Scalar _tmp173 = _tmp140 * _tmp78;
  const Scalar _tmp174 = _tmp140 * _tmp56;
  const Scalar _tmp175 = -_tmp140 * _tmp39 + _tmp172 * _tmp72 - _tmp173 * _tmp38 + _tmp174 * _tmp73;
  const Scalar _tmp176 = _tmp55 * _tmp57;
  const Scalar _tmp177 = _tmp139 * _tmp176;
  const Scalar _tmp178 = _tmp140 * _tmp177;
  const Scalar _tmp179 = _tmp140 * _tmp160;
  const Scalar _tmp180 = -_tmp175 * _tmp84 + _tmp178 * _tmp72 - _tmp179 * _tmp38;
  const Scalar _tmp181 = std::pow(_tmp85, Scalar(-2));
  const Scalar _tmp182 = _tmp180 * _tmp181;
  const Scalar _tmp183 = _tmp108 * _tmp161;
  const Scalar _tmp184 = _tmp56 * _tmp69;
  const Scalar _tmp185 = _tmp140 * _tmp17 - _tmp169 * _tmp184 - _tmp170 * _tmp171 + _tmp173;
  const Scalar _tmp186 = -_tmp169 * _tmp177 + _tmp179 - _tmp185 * _tmp84;
  const Scalar _tmp187 = _tmp108 * _tmp86;
  const Scalar _tmp188 =
      _tmp114 * (-_tmp135 + _tmp136 - _tmp137 + _tmp138 * _tmp66 * _tmp69 - _tmp142 * _tmp170 +
                 _tmp144 - _tmp145 * _tmp162 + _tmp182 * _tmp183 - _tmp186 * _tmp187);
  const Scalar _tmp189 = _tmp147 * _tmp85;
  const Scalar _tmp190 = _tmp163 * _tmp189;
  const Scalar _tmp191 = _tmp164 * _tmp180;
  const Scalar _tmp192 = _tmp186 + _tmp188 - _tmp190 + _tmp191;
  const Scalar _tmp193 = _tmp38 * _tmp86;
  const Scalar _tmp194 = _tmp166 * _tmp182;
  const Scalar _tmp195 = _tmp166 * _tmp86;
  const Scalar _tmp196 = _tmp195 * _tmp38;
  const Scalar _tmp197 = -_tmp140 * _tmp196 - _tmp140 - _tmp167 * _tmp192 + _tmp194 * _tmp72;
  const Scalar _tmp198 = _tmp156 * fh1;
  const Scalar _tmp199 = _tmp103 * _tmp198;
  const Scalar _tmp200 = _tmp157 * (-_tmp149 * _tmp91 + _tmp154);
  const Scalar _tmp201 = _tmp200 * _tmp99;
  const Scalar _tmp202 = _tmp198 * (_tmp119 * _tmp168 + _tmp196 + Scalar(1.0));
  const Scalar _tmp203 = _tmp182 * _tmp88;
  const Scalar _tmp204 = _tmp0 * _tmp103;
  const Scalar _tmp205 = _tmp118 * _tmp150;
  const Scalar _tmp206 = _tmp116 * _tmp182;
  const Scalar _tmp207 = Scalar(1.0) * _tmp172;
  const Scalar _tmp208 = _tmp106 * _tmp174;
  const Scalar _tmp209 = Scalar(1.0) * _tmp178 - _tmp207 * _tmp84 - _tmp208 * _tmp84;
  const Scalar _tmp210 = _tmp108 * _tmp110;
  const Scalar _tmp211 = _tmp114 * (-_tmp106 * _tmp138 - _tmp111 * _tmp145 + Scalar(1.0) * _tmp143 +
                                    _tmp182 * _tmp210 - _tmp187 * _tmp209);
  const Scalar _tmp212 = _tmp113 * _tmp180;
  const Scalar _tmp213 = _tmp112 * _tmp212;
  const Scalar _tmp214 = _tmp112 * _tmp189;
  const Scalar _tmp215 = _tmp209 + _tmp211 + _tmp213 - _tmp214;
  const Scalar _tmp216 = -_tmp120 * _tmp140 - _tmp167 * _tmp215 + _tmp206 * _tmp72;
  const Scalar _tmp217 = _tmp103 * _tmp127;
  const Scalar _tmp218 = _tmp92 * _tmp99;
  const Scalar _tmp219 = _tmp202 * _tmp99;
  const Scalar _tmp220 =
      -_tmp128 * _tmp98 + _tmp129 -
      _tmp158 * (_tmp130 * _tmp148 - _tmp147 * _tmp153 + _tmp149 * _tmp152 - _tmp154 * _tmp155) -
      _tmp199 * (_tmp119 * _tmp197 - _tmp151 * _tmp168 + _tmp192 * _tmp193 - _tmp194 * _tmp38) -
      _tmp200 * _tmp98 + _tmp201 - _tmp202 * _tmp98 -
      _tmp204 * (_tmp152 * _tmp89 - _tmp155 * _tmp90 - _tmp203 * _tmp38 + _tmp203 * _tmp91) -
      _tmp217 * (_tmp119 * _tmp216 - _tmp141 * _tmp205 + _tmp193 * _tmp215 - _tmp206 * _tmp38) +
      _tmp218 + _tmp219 - _tmp92 * _tmp98;
  const Scalar _tmp221 = fh1 * (_tmp40 + _tmp8);
  const Scalar _tmp222 = _tmp123 * fv1 + _tmp156 * _tmp221 + Scalar(40.024799999999999) * _tmp21;
  const Scalar _tmp223 = _tmp34 + _tmp81;
  const Scalar _tmp224 = _tmp223 * _tmp84;
  const Scalar _tmp225 = Scalar(1.0) / (-_tmp224 - _tmp26 + _tmp83);
  const Scalar _tmp226 = Scalar(1.0) * _tmp225;
  const Scalar _tmp227 = _tmp223 * _tmp82;
  const Scalar _tmp228 = _tmp223 * _tmp225;
  const Scalar _tmp229 = _tmp79 * _tmp86;
  const Scalar _tmp230 = -_tmp109 + _tmp115 * _tmp228 - _tmp116 * _tmp229;
  const Scalar _tmp231 = Scalar(1.0) * _tmp82;
  const Scalar _tmp232 = Scalar(1.0) * _tmp127;
  const Scalar _tmp233 = -_tmp121 * fv1 - _tmp126 * _tmp221 - Scalar(40.024799999999999) * _tmp31;
  const Scalar _tmp234 = _tmp226 * _tmp84;
  const Scalar _tmp235 = _tmp224 * _tmp226 + Scalar(1.0);
  const Scalar _tmp236 = _tmp114 * _tmp226;
  const Scalar _tmp237 = -_tmp149 * _tmp79 + _tmp223 * _tmp236;
  const Scalar _tmp238 = Scalar(1.0) * _tmp157;
  const Scalar _tmp239 = _tmp159 + _tmp165 * _tmp228 - _tmp166 * _tmp229;
  const Scalar _tmp240 = Scalar(1.0) * _tmp198;
  const Scalar _tmp241 = _tmp225 * _tmp87;
  const Scalar _tmp242 = -_tmp223 * _tmp241 - _tmp229 * _tmp88 + _tmp81;
  const Scalar _tmp243 =
      Scalar(1.0) * _tmp0 * (-_tmp226 * _tmp87 - _tmp231 * _tmp242 + Scalar(1.0)) +
      Scalar(1.0) * _tmp222 * (_tmp226 * _tmp227 - _tmp226) +
      _tmp232 * (_tmp115 * _tmp226 - _tmp230 * _tmp231) +
      Scalar(1.0) * _tmp233 * (-_tmp231 * _tmp235 + _tmp234) +
      _tmp238 * (-_tmp231 * _tmp237 + _tmp236) + _tmp240 * (_tmp165 * _tmp226 - _tmp231 * _tmp239);
  const Scalar _tmp244 =
      -_tmp103 * _tmp128 - _tmp103 * _tmp200 - _tmp103 * _tmp202 - _tmp103 * _tmp92;
  const Scalar _tmp245 = Scalar(1.0) / (_tmp244);
  const Scalar _tmp246 = std::asinh(_tmp243 * _tmp245);
  const Scalar _tmp247 = Scalar(1.0) * _tmp246;
  const Scalar _tmp248 = Scalar(0.71007031138673404) * _tmp245;
  const Scalar _tmp249 = -_tmp64 + p_a(1, 0);
  const Scalar _tmp250 = -_tmp60 + p_a(0, 0);
  const Scalar _tmp251 =
      std::sqrt(Scalar(std::pow(_tmp249, Scalar(2)) + std::pow(_tmp250, Scalar(2))));
  const Scalar _tmp252 = Scalar(1.4083112389913199) * _tmp246;
  const Scalar _tmp253 = -_tmp244 * _tmp252 - _tmp251;
  const Scalar _tmp254 = _tmp248 * _tmp253;
  const Scalar _tmp255 = Scalar(1.4083112389913199) * _tmp248 * p_a(2, 0) +
                         Scalar(1.4083112389913199) * std::cosh(_tmp247) -
                         Scalar(1.4083112389913199) * std::cosh(_tmp254);
  const Scalar _tmp256 = std::pow(_tmp244, Scalar(-2));
  const Scalar _tmp257 = _tmp243 * _tmp256;
  const Scalar _tmp258 = _tmp79 * _tmp88;
  const Scalar _tmp259 = -_tmp175 * _tmp89 + _tmp182 * _tmp258;
  const Scalar _tmp260 = _tmp0 * _tmp231;
  const Scalar _tmp261 = -_tmp117 * _tmp175 + _tmp206 * _tmp79 + _tmp207 + _tmp208 +
                         _tmp211 * _tmp228 + _tmp213 * _tmp228 - _tmp214 * _tmp228 -
                         _tmp215 * _tmp229;
  const Scalar _tmp262 = _tmp189 * _tmp226;
  const Scalar _tmp263 = _tmp212 * _tmp226;
  const Scalar _tmp264 = Scalar(1.0) * _tmp79;
  const Scalar _tmp265 =
      _tmp147 * _tmp264 - _tmp149 * _tmp175 - _tmp223 * _tmp262 + _tmp223 * _tmp263;
  const Scalar _tmp266 = -_tmp175 * _tmp195 + _tmp185 + _tmp188 * _tmp228 - _tmp190 * _tmp228 +
                         _tmp191 * _tmp228 - _tmp192 * _tmp229 + _tmp194 * _tmp79;
  const Scalar _tmp267 =
      -_tmp220 * _tmp257 + _tmp245 * (_tmp232 * (_tmp211 * _tmp226 + _tmp213 * _tmp226 -
                                                 _tmp214 * _tmp226 - _tmp231 * _tmp261) +
                                      _tmp238 * (-_tmp231 * _tmp265 - _tmp262 + _tmp263) +
                                      _tmp240 * (_tmp188 * _tmp226 - _tmp190 * _tmp226 +
                                                 _tmp191 * _tmp226 - _tmp231 * _tmp266) -
                                      _tmp259 * _tmp260);
  const Scalar _tmp268 =
      std::pow(Scalar(std::pow(_tmp243, Scalar(2)) * _tmp256 + 1), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp269 = Scalar(1.0) * _tmp268 * std::sinh(_tmp247);
  const Scalar _tmp270 = Scalar(1.4083112389913199) * _tmp244;
  const Scalar _tmp271 = _tmp268 * _tmp270;
  const Scalar _tmp272 = Scalar(1.0) / (_tmp251);
  const Scalar _tmp273 = Scalar(0.71007031138673404) * _tmp256;
  const Scalar _tmp274 = _tmp220 * _tmp273;
  const Scalar _tmp275 = std::sinh(_tmp254);
  const Scalar _tmp276 = _tmp0 * _tmp88;
  const Scalar _tmp277 = _tmp182 * _tmp276;
  const Scalar _tmp278 = _tmp127 * _tmp69;
  const Scalar _tmp279 = _tmp118 * _tmp127 * _tmp55;
  const Scalar _tmp280 = _tmp141 * _tmp72;
  const Scalar _tmp281 = _tmp149 * _tmp157;
  const Scalar _tmp282 = _tmp281 * _tmp55;
  const Scalar _tmp283 = _tmp0 * _tmp89;
  const Scalar _tmp284 = _tmp283 * _tmp55;
  const Scalar _tmp285 = _tmp198 * _tmp55;
  const Scalar _tmp286 = _tmp283 * _tmp38;
  const Scalar _tmp287 = _tmp140 * _tmp69;
  const Scalar _tmp288 = _tmp154 * _tmp157;
  const Scalar _tmp289 = _tmp106 * _tmp157;
  const Scalar _tmp290 = _tmp198 * _tmp69;
  const Scalar _tmp291 = -_tmp141 * _tmp168 * _tmp285 - _tmp141 * _tmp279 + _tmp148 * _tmp289 +
                         _tmp197 * _tmp290 + _tmp216 * _tmp278 + _tmp277 * _tmp73 +
                         _tmp280 * _tmp282 + _tmp280 * _tmp284 - _tmp286 * _tmp287 -
                         _tmp287 * _tmp288;
  const Scalar _tmp292 =
      _tmp118 * _tmp278 + _tmp168 * _tmp290 - _tmp281 * _tmp73 - _tmp283 * _tmp73;
  const Scalar _tmp293 = _tmp0 * _tmp82;
  const Scalar _tmp294 = _tmp198 * _tmp82;
  const Scalar _tmp295 = _tmp157 * _tmp82;
  const Scalar _tmp296 = _tmp127 * _tmp82;
  const Scalar _tmp297 = _tmp222 * _tmp226;
  const Scalar _tmp298 = -_tmp227 * _tmp297 + _tmp230 * _tmp296 + _tmp233 * _tmp235 * _tmp82 +
                         _tmp237 * _tmp295 + _tmp239 * _tmp294 + _tmp242 * _tmp293;
  const Scalar _tmp299 = Scalar(1.0) / (_tmp292);
  const Scalar _tmp300 = std::asinh(_tmp298 * _tmp299);
  const Scalar _tmp301 = Scalar(1.4083112389913199) * _tmp300;
  const Scalar _tmp302 =
      -_tmp292 * _tmp301 - std::sqrt(Scalar(std::pow(Scalar(-_tmp45 + p_b(0, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp51 + p_b(1, 0)), Scalar(2))));
  const Scalar _tmp303 = Scalar(0.71007031138673404) * _tmp299;
  const Scalar _tmp304 = _tmp302 * _tmp303;
  const Scalar _tmp305 = Scalar(1.0) * _tmp300;
  const Scalar _tmp306 = Scalar(1.4083112389913199) * _tmp303 * p_b(2, 0) -
                         Scalar(1.4083112389913199) * std::cosh(_tmp304) +
                         Scalar(1.4083112389913199) * std::cosh(_tmp305);
  const Scalar _tmp307 = std::sinh(_tmp304);
  const Scalar _tmp308 = std::pow(_tmp292, Scalar(-2));
  const Scalar _tmp309 = _tmp298 * _tmp308;
  const Scalar _tmp310 = -_tmp291 * _tmp309 + _tmp299 * (_tmp259 * _tmp293 + _tmp261 * _tmp296 +
                                                         _tmp265 * _tmp295 + _tmp266 * _tmp294);
  const Scalar _tmp311 =
      std::pow(Scalar(std::pow(_tmp298, Scalar(2)) * _tmp308 + 1), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp312 = Scalar(1.4083112389913199) * _tmp292;
  const Scalar _tmp313 = _tmp311 * _tmp312;
  const Scalar _tmp314 = Scalar(0.71007031138673404) * _tmp308;
  const Scalar _tmp315 = _tmp302 * _tmp314;
  const Scalar _tmp316 = Scalar(1.0) * _tmp311 * std::sinh(_tmp305);
  const Scalar _tmp317 = _tmp314 * p_b(2, 0);
  const Scalar _tmp318 = _tmp127 * _tmp86;
  const Scalar _tmp319 = _tmp198 * _tmp86;
  const Scalar _tmp320 = -_tmp127 * _tmp206 - _tmp147 * _tmp238 + _tmp192 * _tmp319 -
                         _tmp194 * _tmp198 + _tmp215 * _tmp318 - _tmp277;
  const Scalar _tmp321 = _tmp198 * _tmp225;
  const Scalar _tmp322 = _tmp127 * _tmp225;
  const Scalar _tmp323 = _tmp0 * _tmp241 - _tmp115 * _tmp322 - _tmp157 * _tmp236 -
                         _tmp165 * _tmp321 - _tmp233 * _tmp234 + _tmp297;
  const Scalar _tmp324 = _tmp117 * _tmp127 + _tmp195 * _tmp198 + _tmp281 + _tmp283;
  const Scalar _tmp325 = Scalar(1.0) / (_tmp324);
  const Scalar _tmp326 = std::asinh(_tmp323 * _tmp325);
  const Scalar _tmp327 = Scalar(1.0) * _tmp326;
  const Scalar _tmp328 = Scalar(0.71007031138673404) * _tmp325;
  const Scalar _tmp329 = Scalar(1.4083112389913199) * _tmp324;
  const Scalar _tmp330 =
      -_tmp326 * _tmp329 - std::sqrt(Scalar(std::pow(Scalar(-_tmp27 + p_d(0, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp35 + p_d(1, 0)), Scalar(2))));
  const Scalar _tmp331 = _tmp328 * _tmp330;
  const Scalar _tmp332 = Scalar(1.4083112389913199) * _tmp328 * p_d(2, 0) +
                         Scalar(1.4083112389913199) * std::cosh(_tmp327) -
                         Scalar(1.4083112389913199) * std::cosh(_tmp331);
  const Scalar _tmp333 = std::pow(_tmp324, Scalar(-2));
  const Scalar _tmp334 = Scalar(0.71007031138673404) * _tmp333;
  const Scalar _tmp335 = _tmp320 * _tmp334;
  const Scalar _tmp336 =
      std::pow(Scalar(std::pow(_tmp323, Scalar(2)) * _tmp333 + 1), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp337 = _tmp323 * _tmp333;
  const Scalar _tmp338 =
      _tmp336 *
      (-_tmp320 * _tmp337 +
       _tmp325 * (_tmp157 * _tmp262 - _tmp157 * _tmp263 - _tmp188 * _tmp321 + _tmp190 * _tmp321 -
                  _tmp191 * _tmp321 - _tmp211 * _tmp322 - _tmp213 * _tmp322 + _tmp214 * _tmp322));
  const Scalar _tmp339 = Scalar(1.4083112389913199) * _tmp326;
  const Scalar _tmp340 = std::sinh(_tmp331);
  const Scalar _tmp341 = Scalar(1.0) * std::sinh(_tmp327);
  const Scalar _tmp342 = _tmp139 * _tmp62;
  const Scalar _tmp343 = _tmp342 * _tmp72;
  const Scalar _tmp344 = _tmp62 * _tmp78;
  const Scalar _tmp345 = -_tmp171 * _tmp343 + _tmp344 * _tmp38 + _tmp38 * _tmp75 - _tmp75 * _tmp91;
  const Scalar _tmp346 = _tmp160 * _tmp62;
  const Scalar _tmp347 = -_tmp176 * _tmp343 - _tmp345 * _tmp84 + _tmp346 * _tmp38;
  const Scalar _tmp348 = _tmp113 * _tmp347;
  const Scalar _tmp349 = _tmp163 * _tmp348;
  const Scalar _tmp350 = _tmp140 * _tmp184 + _tmp172 - _tmp344 - _tmp75;
  const Scalar _tmp351 = _tmp178 - _tmp346 - _tmp350 * _tmp84;
  const Scalar _tmp352 = _tmp181 * _tmp347;
  const Scalar _tmp353 = _tmp103 * (-_tmp101 - _tmp132 * _tmp134 + _tmp133 * _tmp95);
  const Scalar _tmp354 = _tmp69 * (-_tmp136 * _tmp67 + _tmp353 * _tmp55);
  const Scalar _tmp355 = _tmp107 * _tmp62;
  const Scalar _tmp356 = -_tmp136 * _tmp71 - _tmp142 * _tmp343 + _tmp353 * _tmp38 -
                         _tmp354 * _tmp72 + _tmp355 * _tmp38;
  const Scalar _tmp357 =
      _tmp114 * (_tmp136 * _tmp66 + _tmp143 - _tmp162 * _tmp356 + _tmp183 * _tmp352 -
                 _tmp187 * _tmp351 - _tmp353 + _tmp354 * _tmp66 - _tmp355);
  const Scalar _tmp358 = _tmp146 * _tmp356;
  const Scalar _tmp359 = _tmp358 * _tmp85;
  const Scalar _tmp360 = _tmp163 * _tmp359;
  const Scalar _tmp361 = _tmp349 + _tmp351 + _tmp357 - _tmp360;
  const Scalar _tmp362 = _tmp168 * _tmp342;
  const Scalar _tmp363 = _tmp352 * _tmp72;
  const Scalar _tmp364 = _tmp166 * _tmp363 - _tmp167 * _tmp361 + _tmp196 * _tmp62 + _tmp62;
  const Scalar _tmp365 = _tmp352 * _tmp38;
  const Scalar _tmp366 = _tmp150 * _tmp343;
  const Scalar _tmp367 = _tmp119 * _tmp62;
  const Scalar _tmp368 = Scalar(1.0) * _tmp342;
  const Scalar _tmp369 = _tmp171 * _tmp368;
  const Scalar _tmp370 = _tmp130 * _tmp75;
  const Scalar _tmp371 = -_tmp176 * _tmp368 + _tmp369 * _tmp84 + _tmp370 * _tmp84;
  const Scalar _tmp372 = _tmp114 * (-_tmp111 * _tmp356 - _tmp142 * _tmp368 - _tmp187 * _tmp371 +
                                    _tmp210 * _tmp352 - Scalar(1.0) * _tmp354);
  const Scalar _tmp373 = _tmp112 * _tmp348;
  const Scalar _tmp374 = _tmp112 * _tmp359;
  const Scalar _tmp375 = _tmp371 + _tmp372 + _tmp373 - _tmp374;
  const Scalar _tmp376 = _tmp116 * _tmp363 + _tmp120 * _tmp62 - _tmp167 * _tmp375;
  const Scalar _tmp377 = _tmp358 * _tmp72;
  const Scalar _tmp378 =
      _tmp129 * _tmp66 -
      _tmp158 * (_tmp130 * _tmp377 - _tmp149 * _tmp366 - _tmp153 * _tmp358 + _tmp154 * _tmp367) -
      _tmp199 * (_tmp119 * _tmp364 + _tmp150 * _tmp362 - _tmp166 * _tmp365 + _tmp193 * _tmp361) +
      _tmp201 * _tmp66 -
      _tmp204 *
          (_tmp352 * _tmp88 * _tmp91 - _tmp365 * _tmp88 - _tmp366 * _tmp89 + _tmp367 * _tmp90) -
      _tmp217 * (-_tmp116 * _tmp365 + _tmp119 * _tmp376 + _tmp193 * _tmp375 + _tmp205 * _tmp342) +
      _tmp218 * _tmp66 + _tmp219 * _tmp66;
  const Scalar _tmp379 = _tmp352 * _tmp79;
  const Scalar _tmp380 = _tmp166 * _tmp379 - _tmp195 * _tmp345 + _tmp228 * _tmp349 +
                         _tmp228 * _tmp357 - _tmp228 * _tmp360 - _tmp229 * _tmp361 + _tmp350;
  const Scalar _tmp381 = _tmp226 * _tmp348;
  const Scalar _tmp382 = _tmp226 * _tmp359;
  const Scalar _tmp383 =
      -_tmp149 * _tmp345 + _tmp223 * _tmp381 - _tmp223 * _tmp382 + _tmp264 * _tmp358;
  const Scalar _tmp384 = _tmp258 * _tmp352 - _tmp345 * _tmp89;
  const Scalar _tmp385 = _tmp116 * _tmp379 - _tmp117 * _tmp345 + _tmp228 * _tmp372 +
                         _tmp228 * _tmp373 - _tmp228 * _tmp374 - _tmp229 * _tmp375 - _tmp369 -
                         _tmp370;
  const Scalar _tmp386 = _tmp245 * (_tmp232 * (_tmp226 * _tmp372 + _tmp226 * _tmp373 -
                                               _tmp226 * _tmp374 - _tmp231 * _tmp385) +
                                    _tmp238 * (-_tmp231 * _tmp383 + _tmp381 - _tmp382) +
                                    _tmp240 * (_tmp226 * _tmp349 + _tmp226 * _tmp357 -
                                               _tmp226 * _tmp360 - _tmp231 * _tmp380) -
                                    _tmp260 * _tmp384) -
                         _tmp257 * _tmp378;
  const Scalar _tmp387 = _tmp273 * _tmp378;
  const Scalar _tmp388 = _tmp276 * _tmp352;
  const Scalar _tmp389 = _tmp62 * _tmp69;
  const Scalar _tmp390 = _tmp278 * _tmp376 + _tmp279 * _tmp342 - _tmp282 * _tmp343 -
                         _tmp284 * _tmp343 + _tmp285 * _tmp362 + _tmp286 * _tmp389 +
                         _tmp288 * _tmp389 + _tmp289 * _tmp377 + _tmp290 * _tmp364 +
                         _tmp388 * _tmp73;
  const Scalar _tmp391 =
      _tmp299 * (_tmp293 * _tmp384 + _tmp294 * _tmp380 + _tmp295 * _tmp383 + _tmp296 * _tmp385) -
      _tmp309 * _tmp390;
  const Scalar _tmp392 = -_tmp116 * _tmp127 * _tmp352 - _tmp166 * _tmp198 * _tmp352 -
                         _tmp238 * _tmp358 + _tmp318 * _tmp375 + _tmp319 * _tmp361 - _tmp388;
  const Scalar _tmp393 = _tmp334 * _tmp392;
  const Scalar _tmp394 =
      _tmp336 *
      (_tmp325 * (-_tmp157 * _tmp381 + _tmp157 * _tmp382 - _tmp321 * _tmp349 - _tmp321 * _tmp357 +
                  _tmp321 * _tmp360 - _tmp322 * _tmp372 - _tmp322 * _tmp373 + _tmp322 * _tmp374) -
       _tmp337 * _tmp392);

  // Output terms (1)
  Eigen::Matrix<Scalar, 4, 3> _res;

  _res(0, 0) = 0;
  _res(1, 0) =
      -_tmp220 * _tmp255 -
      _tmp270 * (_tmp267 * _tmp269 - _tmp274 * p_a(2, 0) -
                 _tmp275 * (_tmp248 * (-_tmp220 * _tmp252 - _tmp250 * _tmp272 - _tmp267 * _tmp271) -
                            _tmp253 * _tmp274));
  _res(2, 0) =
      -_tmp291 * _tmp306 -
      _tmp312 *
          (-_tmp291 * _tmp317 -
           _tmp307 * (-_tmp291 * _tmp315 + _tmp303 * (-_tmp291 * _tmp301 - _tmp310 * _tmp313)) +
           _tmp310 * _tmp316);
  _res(3, 0) =
      -_tmp320 * _tmp332 -
      _tmp329 *
          (-_tmp335 * p_d(2, 0) + _tmp338 * _tmp341 -
           _tmp340 * (_tmp328 * (-_tmp320 * _tmp339 - _tmp329 * _tmp338) - _tmp330 * _tmp335));
  _res(0, 1) = 0;
  _res(1, 1) =
      -_tmp255 * _tmp378 -
      _tmp270 * (_tmp269 * _tmp386 -
                 _tmp275 * (_tmp248 * (-_tmp249 * _tmp272 - _tmp252 * _tmp378 - _tmp271 * _tmp386) -
                            _tmp253 * _tmp387) -
                 _tmp387 * p_a(2, 0));
  _res(2, 1) =
      -_tmp306 * _tmp390 -
      _tmp312 *
          (-_tmp307 * (_tmp303 * (-_tmp301 * _tmp390 - _tmp313 * _tmp391) - _tmp315 * _tmp390) +
           _tmp316 * _tmp391 - _tmp317 * _tmp390);
  _res(3, 1) =
      -_tmp329 *
          (-_tmp340 * (_tmp328 * (-_tmp329 * _tmp394 - _tmp339 * _tmp392) - _tmp330 * _tmp393) +
           _tmp341 * _tmp394 - _tmp393 * p_d(2, 0)) -
      _tmp332 * _tmp392;
  _res(0, 2) = 0;
  _res(1, 2) = Scalar(-1.0);
  _res(2, 2) = 0;
  _res(3, 2) = 0;

  return _res;
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym

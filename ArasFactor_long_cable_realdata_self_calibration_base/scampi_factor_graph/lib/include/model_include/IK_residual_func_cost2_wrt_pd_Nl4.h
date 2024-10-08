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
 * Symbolic function: IK_residual_func_cost2_wrt_pd_Nl4
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
Eigen::Matrix<Scalar, 4, 3> IkResidualFuncCost2WrtPdNl4(
    const Scalar fh1, const Scalar fv1, const sym::Rot3<Scalar>& DeltaRot,
    const Eigen::Matrix<Scalar, 3, 1>& position_vector, const Eigen::Matrix<Scalar, 4, 1>& encoder,
    const Eigen::Matrix<Scalar, 3, 1>& p_a, const Eigen::Matrix<Scalar, 3, 1>& p_b,
    const Eigen::Matrix<Scalar, 3, 1>& p_c, const Eigen::Matrix<Scalar, 3, 1>& p_d,
    const sym::Rot3<Scalar>& Rot_init, const Scalar epsilon) {
  // Total ops: 1212

  // Unused inputs
  (void)encoder;
  (void)epsilon;

  // Input arrays
  const Eigen::Matrix<Scalar, 4, 1>& _DeltaRot = DeltaRot.Data();
  const Eigen::Matrix<Scalar, 4, 1>& _Rot_init = Rot_init.Data();

  // Intermediate terms (400)
  const Scalar _tmp0 = _DeltaRot[0] * _Rot_init[3] - _DeltaRot[1] * _Rot_init[2] +
                       _DeltaRot[2] * _Rot_init[1] + _DeltaRot[3] * _Rot_init[0];
  const Scalar _tmp1 = _DeltaRot[0] * _Rot_init[2] + _DeltaRot[1] * _Rot_init[3] -
                       _DeltaRot[2] * _Rot_init[0] + _DeltaRot[3] * _Rot_init[1];
  const Scalar _tmp2 = 2 * _tmp1;
  const Scalar _tmp3 = _tmp0 * _tmp2;
  const Scalar _tmp4 = -_DeltaRot[0] * _Rot_init[1] + _DeltaRot[1] * _Rot_init[0] +
                       _DeltaRot[2] * _Rot_init[3] + _DeltaRot[3] * _Rot_init[2];
  const Scalar _tmp5 = -2 * _DeltaRot[0] * _Rot_init[0] - 2 * _DeltaRot[1] * _Rot_init[1] -
                       2 * _DeltaRot[2] * _Rot_init[2] + 2 * _DeltaRot[3] * _Rot_init[3];
  const Scalar _tmp6 = _tmp4 * _tmp5;
  const Scalar _tmp7 = Scalar(0.20999999999999999) * _tmp3 - Scalar(0.20999999999999999) * _tmp6;
  const Scalar _tmp8 = 2 * _tmp0 * _tmp4;
  const Scalar _tmp9 = _tmp1 * _tmp5;
  const Scalar _tmp10 = _tmp8 + _tmp9;
  const Scalar _tmp11 = -Scalar(0.010999999999999999) * _tmp10;
  const Scalar _tmp12 = -2 * std::pow(_tmp1, Scalar(2));
  const Scalar _tmp13 = 1 - 2 * std::pow(_tmp4, Scalar(2));
  const Scalar _tmp14 = Scalar(0.20999999999999999) * _tmp12 + Scalar(0.20999999999999999) * _tmp13;
  const Scalar _tmp15 = _tmp11 - _tmp14;
  const Scalar _tmp16 = _tmp15 + _tmp7;
  const Scalar _tmp17 = _tmp16 + position_vector(0, 0);
  const Scalar _tmp18 = -_tmp17 + p_d(0, 0);
  const Scalar _tmp19 = _tmp2 * _tmp4;
  const Scalar _tmp20 = _tmp0 * _tmp5;
  const Scalar _tmp21 = Scalar(0.20999999999999999) * _tmp19 + Scalar(0.20999999999999999) * _tmp20;
  const Scalar _tmp22 = -2 * std::pow(_tmp0, Scalar(2));
  const Scalar _tmp23 = -Scalar(0.010999999999999999) * _tmp12 -
                        Scalar(0.010999999999999999) * _tmp22 + Scalar(-0.010999999999999999);
  const Scalar _tmp24 = Scalar(0.20999999999999999) * _tmp8 - Scalar(0.20999999999999999) * _tmp9;
  const Scalar _tmp25 = _tmp23 - _tmp24;
  const Scalar _tmp26 = _tmp21 + _tmp25;
  const Scalar _tmp27 = -_tmp26 + p_d(2, 0) - position_vector(2, 0);
  const Scalar _tmp28 = Scalar(0.20999999999999999) * _tmp13 + Scalar(0.20999999999999999) * _tmp22;
  const Scalar _tmp29 = Scalar(0.20999999999999999) * _tmp3 + Scalar(0.20999999999999999) * _tmp6;
  const Scalar _tmp30 = _tmp19 - _tmp20;
  const Scalar _tmp31 = -Scalar(0.010999999999999999) * _tmp30;
  const Scalar _tmp32 = -_tmp29 + _tmp31;
  const Scalar _tmp33 = _tmp28 + _tmp32;
  const Scalar _tmp34 = _tmp33 + position_vector(1, 0);
  const Scalar _tmp35 = -_tmp34 + p_d(1, 0);
  const Scalar _tmp36 = std::pow(_tmp18, Scalar(2)) + std::pow(_tmp35, Scalar(2));
  const Scalar _tmp37 =
      std::pow(Scalar(std::pow(_tmp27, Scalar(2)) + _tmp36), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp38 = -_tmp7;
  const Scalar _tmp39 = _tmp15 + _tmp38;
  const Scalar _tmp40 = _tmp39 - p_a(0, 0) + position_vector(0, 0);
  const Scalar _tmp41 = -_tmp28;
  const Scalar _tmp42 = _tmp32 + _tmp41;
  const Scalar _tmp43 = _tmp42 - p_a(1, 0) + position_vector(1, 0);
  const Scalar _tmp44 = std::pow(Scalar(std::pow(_tmp40, Scalar(2)) + std::pow(_tmp43, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp45 = _tmp43 * _tmp44;
  const Scalar _tmp46 = -_tmp21;
  const Scalar _tmp47 = fh1 * (_tmp25 + _tmp46);
  const Scalar _tmp48 = -Scalar(40.024799999999999) * _tmp30 - _tmp42 * fv1 - _tmp45 * _tmp47;
  const Scalar _tmp49 = Scalar(1.0) * _tmp16;
  const Scalar _tmp50 = Scalar(1.0) * _tmp33;
  const Scalar _tmp51 = -_tmp50;
  const Scalar _tmp52 = _tmp29 + _tmp31;
  const Scalar _tmp53 = _tmp28 + _tmp52;
  const Scalar _tmp54 = _tmp51 + _tmp53;
  const Scalar _tmp55 = _tmp11 + _tmp14;
  const Scalar _tmp56 = _tmp38 + _tmp55;
  const Scalar _tmp57 = _tmp41 + _tmp52;
  const Scalar _tmp58 = Scalar(1.0) / (_tmp51 + _tmp57);
  const Scalar _tmp59 = _tmp58 * (_tmp49 - _tmp56);
  const Scalar _tmp60 = _tmp54 * _tmp59;
  const Scalar _tmp61 = _tmp55 + _tmp7;
  const Scalar _tmp62 = Scalar(1.0) / (_tmp49 - _tmp60 - _tmp61);
  const Scalar _tmp63 = Scalar(1.0) * _tmp62;
  const Scalar _tmp64 = _tmp60 * _tmp63 + Scalar(1.0);
  const Scalar _tmp65 = Scalar(1.0) * _tmp58;
  const Scalar _tmp66 = _tmp59 * _tmp63;
  const Scalar _tmp67 = _tmp17 - p_d(0, 0);
  const Scalar _tmp68 = Scalar(1.0) / (_tmp67);
  const Scalar _tmp69 = _tmp34 - p_d(1, 0);
  const Scalar _tmp70 = std::pow(_tmp67, Scalar(2));
  const Scalar _tmp71 = std::pow(_tmp69, Scalar(2));
  const Scalar _tmp72 = _tmp70 + _tmp71;
  const Scalar _tmp73 = std::sqrt(_tmp72);
  const Scalar _tmp74 = Scalar(1.0) / (_tmp73);
  const Scalar _tmp75 = _tmp16 * _tmp74;
  const Scalar _tmp76 = _tmp33 * _tmp74;
  const Scalar _tmp77 = -_tmp67 * _tmp76 + _tmp69 * _tmp75;
  const Scalar _tmp78 = _tmp73 * _tmp77;
  const Scalar _tmp79 = _tmp68 * _tmp78;
  const Scalar _tmp80 = _tmp61 + position_vector(0, 0);
  const Scalar _tmp81 = _tmp80 - p_c(0, 0);
  const Scalar _tmp82 = _tmp53 + position_vector(1, 0);
  const Scalar _tmp83 = _tmp82 - p_c(1, 0);
  const Scalar _tmp84 = std::pow(Scalar(std::pow(_tmp81, Scalar(2)) + std::pow(_tmp83, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp85 = _tmp81 * _tmp84;
  const Scalar _tmp86 = _tmp83 * _tmp84;
  const Scalar _tmp87 = _tmp56 + position_vector(0, 0);
  const Scalar _tmp88 = _tmp87 - p_b(0, 0);
  const Scalar _tmp89 = std::pow(_tmp88, Scalar(2));
  const Scalar _tmp90 = _tmp57 + position_vector(1, 0);
  const Scalar _tmp91 = _tmp90 - p_b(1, 0);
  const Scalar _tmp92 = _tmp89 + std::pow(_tmp91, Scalar(2));
  const Scalar _tmp93 = std::pow(_tmp92, Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp94 = _tmp91 * _tmp93;
  const Scalar _tmp95 = _tmp88 * _tmp93;
  const Scalar _tmp96 = -_tmp56 * _tmp94 + _tmp57 * _tmp95 + _tmp79 * _tmp95;
  const Scalar _tmp97 = _tmp68 * _tmp69;
  const Scalar _tmp98 = _tmp85 * _tmp97;
  const Scalar _tmp99 = -_tmp86 + _tmp98;
  const Scalar _tmp100 = _tmp95 * _tmp97;
  const Scalar _tmp101 = _tmp100 - _tmp94;
  const Scalar _tmp102 = Scalar(1.0) / (_tmp101);
  const Scalar _tmp103 = _tmp102 * _tmp99;
  const Scalar _tmp104 = -_tmp103 * _tmp96 + _tmp53 * _tmp85 - _tmp61 * _tmp86 + _tmp79 * _tmp85;
  const Scalar _tmp105 = _tmp23 + _tmp24;
  const Scalar _tmp106 = _tmp105 + _tmp46;
  const Scalar _tmp107 = _tmp26 * _tmp68;
  const Scalar _tmp108 = _tmp107 * _tmp69;
  const Scalar _tmp109 = _tmp106 * _tmp94 - _tmp108 * _tmp95;
  const Scalar _tmp110 = _tmp102 * _tmp109;
  const Scalar _tmp111 = Scalar(1.0) * _tmp110;
  const Scalar _tmp112 = _tmp26 * _tmp95;
  const Scalar _tmp113 = -_tmp106 * _tmp95 + _tmp112;
  const Scalar _tmp114 = Scalar(1.0) * _tmp102;
  const Scalar _tmp115 = _tmp111 * _tmp59 - _tmp113 * _tmp114;
  const Scalar _tmp116 = _tmp105 + _tmp21;
  const Scalar _tmp117 = -_tmp108 * _tmp85 - _tmp110 * _tmp99 + _tmp116 * _tmp86;
  const Scalar _tmp118 = -_tmp103 * _tmp113 - _tmp116 * _tmp85 - _tmp117 * _tmp59 + _tmp26 * _tmp85;
  const Scalar _tmp119 = Scalar(1.0) / (_tmp118);
  const Scalar _tmp120 = _tmp115 * _tmp119;
  const Scalar _tmp121 = -_tmp104 * _tmp120 - _tmp114 * _tmp96;
  const Scalar _tmp122 = Scalar(1.0) / (_tmp104);
  const Scalar _tmp123 = _tmp118 * _tmp122;
  const Scalar _tmp124 = _tmp121 * _tmp123;
  const Scalar _tmp125 = _tmp115 + _tmp124;
  const Scalar _tmp126 = _tmp117 * _tmp119;
  const Scalar _tmp127 = _tmp54 * _tmp62;
  const Scalar _tmp128 = -_tmp111 + _tmp124 * _tmp127 - _tmp125 * _tmp126;
  const Scalar _tmp129 = _tmp45 * fh1;
  const Scalar _tmp130 = Scalar(1.0) * _tmp129;
  const Scalar _tmp131 = _tmp40 * _tmp44;
  const Scalar _tmp132 = Scalar(40.024799999999999) * _tmp10 + _tmp131 * _tmp47 + _tmp39 * fv1;
  const Scalar _tmp133 = _tmp54 * _tmp63;
  const Scalar _tmp134 = _tmp102 * _tmp96;
  const Scalar _tmp135 = _tmp102 * _tmp113;
  const Scalar _tmp136 = _tmp108 + _tmp110 * _tmp97;
  const Scalar _tmp137 = _tmp135 * _tmp97 - _tmp136 * _tmp59 - _tmp26;
  const Scalar _tmp138 = _tmp104 * _tmp119;
  const Scalar _tmp139 = _tmp134 * _tmp97 - _tmp137 * _tmp138 - _tmp79;
  const Scalar _tmp140 = _tmp123 * _tmp139;
  const Scalar _tmp141 = _tmp137 + _tmp140;
  const Scalar _tmp142 = -_tmp126 * _tmp141 + _tmp127 * _tmp140 + _tmp136;
  const Scalar _tmp143 = _tmp131 * fh1;
  const Scalar _tmp144 = Scalar(1.0) * _tmp143;
  const Scalar _tmp145 = _tmp123 * _tmp63;
  const Scalar _tmp146 = Scalar(1.0) * _tmp122;
  const Scalar _tmp147 = -_tmp117 * _tmp146 + _tmp123 * _tmp133;
  const Scalar _tmp148 = fh1 * (-_tmp131 * _tmp42 + _tmp39 * _tmp45);
  const Scalar _tmp149 = Scalar(1.0) * _tmp148;
  const Scalar _tmp150 = _tmp49 + _tmp50 * _tmp59;
  const Scalar _tmp151 = 0;
  const Scalar _tmp152 = _tmp150 * _tmp62;
  const Scalar _tmp153 = -_tmp126 * _tmp151 - _tmp152 * _tmp54 + _tmp51;
  const Scalar _tmp154 = Scalar(333.54000000000002) - fv1;
  const Scalar _tmp155 = Scalar(1.0) * _tmp154;
  const Scalar _tmp156 = _tmp130 * (_tmp124 * _tmp63 - _tmp128 * _tmp65) +
                         Scalar(1.0) * _tmp132 * (_tmp133 * _tmp58 - _tmp63) +
                         _tmp144 * (_tmp140 * _tmp63 - _tmp142 * _tmp65) +
                         _tmp149 * (_tmp145 - _tmp147 * _tmp65) +
                         _tmp155 * (-_tmp150 * _tmp63 - _tmp153 * _tmp65 + Scalar(1.0)) +
                         Scalar(1.0) * _tmp48 * (-_tmp64 * _tmp65 + _tmp66);
  const Scalar _tmp157 = _tmp68 * _tmp73;
  const Scalar _tmp158 = _tmp102 * _tmp95;
  const Scalar _tmp159 = _tmp119 * _tmp99;
  const Scalar _tmp160 = _tmp151 * _tmp159;
  const Scalar _tmp161 = _tmp119 * _tmp151;
  const Scalar _tmp162 = _tmp161 * _tmp85;
  const Scalar _tmp163 = _tmp154 * (-_tmp158 * _tmp160 + _tmp162);
  const Scalar _tmp164 = -_tmp125 * _tmp159 + Scalar(1.0);
  const Scalar _tmp165 = _tmp119 * _tmp85;
  const Scalar _tmp166 = _tmp125 * _tmp165;
  const Scalar _tmp167 = _tmp129 * (_tmp158 * _tmp164 + _tmp166);
  const Scalar _tmp168 = _tmp141 * _tmp165;
  const Scalar _tmp169 = -_tmp141 * _tmp159 - _tmp97;
  const Scalar _tmp170 = _tmp143 * (_tmp158 * _tmp169 + _tmp168 + Scalar(1.0));
  const Scalar _tmp171 = _tmp146 * _tmp85;
  const Scalar _tmp172 = _tmp103 * _tmp95;
  const Scalar _tmp173 = _tmp148 * (-_tmp146 * _tmp172 + _tmp171);
  const Scalar _tmp174 =
      -_tmp157 * _tmp163 - _tmp157 * _tmp167 - _tmp157 * _tmp170 - _tmp157 * _tmp173;
  const Scalar _tmp175 = Scalar(1.0) / (_tmp174);
  const Scalar _tmp176 = std::asinh(_tmp156 * _tmp175);
  const Scalar _tmp177 = Scalar(1.0) * _tmp176;
  const Scalar _tmp178 = Scalar(1.0) * std::cosh(_tmp177);
  const Scalar _tmp179 = std::pow(_tmp174, Scalar(-2));
  const Scalar _tmp180 =
      std::pow(Scalar(std::pow(_tmp156, Scalar(2)) * _tmp179 + 1), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp181 = Scalar(1.0) / (_tmp70);
  const Scalar _tmp182 = _tmp181 * _tmp69;
  const Scalar _tmp183 = _tmp135 * _tmp182;
  const Scalar _tmp184 = _tmp110 * _tmp182;
  const Scalar _tmp185 = std::pow(_tmp101, Scalar(-2));
  const Scalar _tmp186 = _tmp185 * _tmp95;
  const Scalar _tmp187 = _tmp109 * _tmp186;
  const Scalar _tmp188 = _tmp182 * _tmp187;
  const Scalar _tmp189 = _tmp182 * _tmp26;
  const Scalar _tmp190 = _tmp172 * _tmp189 - _tmp184 * _tmp85 + _tmp188 * _tmp99 - _tmp189 * _tmp85;
  const Scalar _tmp191 = _tmp113 * _tmp186;
  const Scalar _tmp192 = _tmp182 * _tmp191;
  const Scalar _tmp193 = -_tmp183 * _tmp85 - _tmp190 * _tmp59 + _tmp192 * _tmp99;
  const Scalar _tmp194 = _tmp122 * _tmp193;
  const Scalar _tmp195 = _tmp74 * _tmp77;
  const Scalar _tmp196 = _tmp134 * _tmp182;
  const Scalar _tmp197 = _tmp186 * _tmp96;
  const Scalar _tmp198 = _tmp182 * _tmp197;
  const Scalar _tmp199 = _tmp181 * _tmp78;
  const Scalar _tmp200 = std::pow(_tmp72, Scalar(Scalar(-3) / Scalar(2)));
  const Scalar _tmp201 = _tmp200 * _tmp33;
  const Scalar _tmp202 = _tmp16 * _tmp200;
  const Scalar _tmp203 = _tmp67 * _tmp69;
  const Scalar _tmp204 = _tmp157 * (-_tmp201 * _tmp70 + _tmp202 * _tmp203 + _tmp76);
  const Scalar _tmp205 = -_tmp195 * _tmp95 + _tmp199 * _tmp95 + _tmp204 * _tmp95;
  const Scalar _tmp206 = -_tmp103 * _tmp205 - _tmp195 * _tmp85 - _tmp196 * _tmp85 +
                         _tmp198 * _tmp99 + _tmp199 * _tmp85 + _tmp204 * _tmp85;
  const Scalar _tmp207 = std::pow(_tmp104, Scalar(-2));
  const Scalar _tmp208 = _tmp206 * _tmp207;
  const Scalar _tmp209 = _tmp118 * _tmp133;
  const Scalar _tmp210 = Scalar(1.0) * _tmp117;
  const Scalar _tmp211 =
      _tmp133 * _tmp194 - _tmp146 * _tmp190 - _tmp208 * _tmp209 + _tmp208 * _tmp210;
  const Scalar _tmp212 = _tmp118 * _tmp63;
  const Scalar _tmp213 = _tmp208 * _tmp212;
  const Scalar _tmp214 = _tmp194 * _tmp63;
  const Scalar _tmp215 = _tmp119 * _tmp190;
  const Scalar _tmp216 = std::pow(_tmp118, Scalar(-2));
  const Scalar _tmp217 = _tmp193 * _tmp216;
  const Scalar _tmp218 = _tmp117 * _tmp217;
  const Scalar _tmp219 = -_tmp151 * _tmp215 + _tmp151 * _tmp218;
  const Scalar _tmp220 = _tmp155 * _tmp58;
  const Scalar _tmp221 = _tmp71 / [&]() {
    const Scalar base = _tmp67;
    return base * base * base;
  }();
  const Scalar _tmp222 = _tmp102 * _tmp97;
  const Scalar _tmp223 = _tmp119 * _tmp137;
  const Scalar _tmp224 = -_tmp102 * _tmp112 * _tmp221 + _tmp184 - _tmp187 * _tmp221 + _tmp189;
  const Scalar _tmp225 = _tmp183 - _tmp191 * _tmp221 - _tmp224 * _tmp59;
  const Scalar _tmp226 = _tmp104 * _tmp217;
  const Scalar _tmp227 =
      _tmp123 * (_tmp137 * _tmp226 - _tmp138 * _tmp225 + _tmp195 + _tmp196 - _tmp197 * _tmp221 -
                 _tmp199 - _tmp204 + _tmp205 * _tmp222 - _tmp206 * _tmp223);
  const Scalar _tmp228 = _tmp118 * _tmp139;
  const Scalar _tmp229 = _tmp208 * _tmp228;
  const Scalar _tmp230 = _tmp139 * _tmp194;
  const Scalar _tmp231 = _tmp225 + _tmp227 - _tmp229 + _tmp230;
  const Scalar _tmp232 = _tmp141 * _tmp216;
  const Scalar _tmp233 = _tmp193 * _tmp232;
  const Scalar _tmp234 = _tmp117 * _tmp233 - _tmp126 * _tmp231 + _tmp127 * _tmp227 -
                         _tmp127 * _tmp229 + _tmp127 * _tmp230 - _tmp141 * _tmp215 + _tmp224;
  const Scalar _tmp235 = _tmp118 * _tmp121;
  const Scalar _tmp236 = _tmp208 * _tmp235;
  const Scalar _tmp237 = _tmp121 * _tmp194;
  const Scalar _tmp238 = _tmp114 * _tmp95;
  const Scalar _tmp239 = _tmp189 * _tmp238;
  const Scalar _tmp240 = Scalar(1.0) * _tmp109;
  const Scalar _tmp241 = _tmp182 * _tmp186;
  const Scalar _tmp242 = _tmp240 * _tmp241;
  const Scalar _tmp243 = Scalar(1.0) * _tmp192 - _tmp239 * _tmp59 - _tmp242 * _tmp59;
  const Scalar _tmp244 = _tmp123 * (-_tmp114 * _tmp205 + _tmp115 * _tmp226 - _tmp120 * _tmp206 -
                                    _tmp138 * _tmp243 + Scalar(1.0) * _tmp198);
  const Scalar _tmp245 = -_tmp236 + _tmp237 + _tmp243 + _tmp244;
  const Scalar _tmp246 = -_tmp125 * _tmp215 + _tmp125 * _tmp218 - _tmp126 * _tmp245 -
                         _tmp127 * _tmp236 + _tmp127 * _tmp237 + _tmp127 * _tmp244 + _tmp239 +
                         _tmp242;
  const Scalar _tmp247 = _tmp181 * _tmp73;
  const Scalar _tmp248 = _tmp173 * _tmp74;
  const Scalar _tmp249 = _tmp125 * _tmp99;
  const Scalar _tmp250 = -_tmp159 * _tmp245 - _tmp166 * _tmp182 + _tmp217 * _tmp249;
  const Scalar _tmp251 = _tmp125 * _tmp85;
  const Scalar _tmp252 = _tmp89 / _tmp92;
  const Scalar _tmp253 = _tmp182 * _tmp185;
  const Scalar _tmp254 = _tmp252 * _tmp253;
  const Scalar _tmp255 = _tmp129 * _tmp157;
  const Scalar _tmp256 = _tmp167 * _tmp74;
  const Scalar _tmp257 = _tmp208 * _tmp99;
  const Scalar _tmp258 = Scalar(1.0) * _tmp85;
  const Scalar _tmp259 = _tmp146 * _tmp99;
  const Scalar _tmp260 = _tmp158 * _tmp182;
  const Scalar _tmp261 = _tmp148 * _tmp157;
  const Scalar _tmp262 = _tmp170 * _tmp74;
  const Scalar _tmp263 = _tmp232 * _tmp99;
  const Scalar _tmp264 = -_tmp159 * _tmp231 - _tmp168 * _tmp182 - _tmp182 + _tmp193 * _tmp263;
  const Scalar _tmp265 = _tmp169 * _tmp252;
  const Scalar _tmp266 = _tmp143 * _tmp157;
  const Scalar _tmp267 = _tmp163 * _tmp74;
  const Scalar _tmp268 = _tmp151 * _tmp172;
  const Scalar _tmp269 = _tmp151 * _tmp85;
  const Scalar _tmp270 = _tmp154 * _tmp157;
  const Scalar _tmp271 =
      -_tmp163 * _tmp247 - _tmp167 * _tmp247 - _tmp170 * _tmp247 - _tmp173 * _tmp247 + _tmp248 -
      _tmp255 * (_tmp158 * _tmp250 - _tmp164 * _tmp254 + _tmp165 * _tmp245 - _tmp217 * _tmp251) +
      _tmp256 -
      _tmp261 * (-_tmp171 * _tmp260 - _tmp208 * _tmp258 + _tmp238 * _tmp257 + _tmp254 * _tmp259) +
      _tmp262 -
      _tmp266 * (_tmp158 * _tmp264 + _tmp165 * _tmp231 - _tmp233 * _tmp85 - _tmp253 * _tmp265) +
      _tmp267 -
      _tmp270 * (_tmp160 * _tmp254 - _tmp162 * _tmp260 + _tmp217 * _tmp268 - _tmp217 * _tmp269);
  const Scalar _tmp272 = _tmp156 * _tmp179;
  const Scalar _tmp273 =
      _tmp180 *
      (_tmp175 *
           (_tmp130 * (-_tmp236 * _tmp63 + _tmp237 * _tmp63 + _tmp244 * _tmp63 - _tmp246 * _tmp65) +
            _tmp144 * (_tmp227 * _tmp63 - _tmp229 * _tmp63 + _tmp230 * _tmp63 - _tmp234 * _tmp65) +
            _tmp149 * (-_tmp211 * _tmp65 - _tmp213 + _tmp214) - _tmp219 * _tmp220) -
       _tmp271 * _tmp272);
  const Scalar _tmp274 = std::sqrt(_tmp36);
  const Scalar _tmp275 = Scalar(1.4083112389913199) * _tmp174;
  const Scalar _tmp276 = -_tmp176 * _tmp275 - _tmp274;
  const Scalar _tmp277 = Scalar(0.71007031138673404) * _tmp179 * _tmp276;
  const Scalar _tmp278 = Scalar(1.4083112389913199) * _tmp271;
  const Scalar _tmp279 = Scalar(1.0) / (_tmp274);
  const Scalar _tmp280 = Scalar(0.71007031138673404) * _tmp175;
  const Scalar _tmp281 = _tmp276 * _tmp280;
  const Scalar _tmp282 = std::cosh(_tmp281);
  const Scalar _tmp283 = -std::sinh(_tmp177) - std::sinh(_tmp281);
  const Scalar _tmp284 = _tmp151 * _tmp154;
  const Scalar _tmp285 = _tmp159 * _tmp284;
  const Scalar _tmp286 = _tmp146 * _tmp148;
  const Scalar _tmp287 = _tmp102 * _tmp143;
  const Scalar _tmp288 = _tmp102 * _tmp129;
  const Scalar _tmp289 =
      -_tmp102 * _tmp285 - _tmp103 * _tmp286 + _tmp164 * _tmp288 + _tmp169 * _tmp287;
  const Scalar _tmp290 = Scalar(1.0) / (_tmp289);
  const Scalar _tmp291 = _tmp154 * _tmp58;
  const Scalar _tmp292 = _tmp129 * _tmp58;
  const Scalar _tmp293 = _tmp143 * _tmp58;
  const Scalar _tmp294 = _tmp148 * _tmp58;
  const Scalar _tmp295 = _tmp132 * _tmp63;
  const Scalar _tmp296 = _tmp128 * _tmp292 + _tmp142 * _tmp293 + _tmp147 * _tmp294 +
                         _tmp153 * _tmp291 - _tmp295 * _tmp54 * _tmp58 + _tmp48 * _tmp58 * _tmp64;
  const Scalar _tmp297 = std::asinh(_tmp290 * _tmp296);
  const Scalar _tmp298 = Scalar(1.4083112389913199) * _tmp297;
  const Scalar _tmp299 =
      -_tmp289 * _tmp298 - std::sqrt(Scalar(std::pow(Scalar(-_tmp87 + p_b(0, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp90 + p_b(1, 0)), Scalar(2))));
  const Scalar _tmp300 = Scalar(0.71007031138673404) * _tmp290;
  const Scalar _tmp301 = _tmp299 * _tmp300;
  const Scalar _tmp302 = Scalar(1.0) * _tmp297;
  const Scalar _tmp303 = -std::sinh(_tmp301) - std::sinh(_tmp302);
  const Scalar _tmp304 = _tmp114 * _tmp148;
  const Scalar _tmp305 = _tmp217 * _tmp284;
  const Scalar _tmp306 = _tmp143 * _tmp169;
  const Scalar _tmp307 = _tmp154 * _tmp161;
  const Scalar _tmp308 = _tmp307 * _tmp85;
  const Scalar _tmp309 = _tmp102 * _tmp182;
  const Scalar _tmp310 = _tmp286 * _tmp85;
  const Scalar _tmp311 = _tmp103 * _tmp305 - _tmp129 * _tmp164 * _tmp241 + _tmp241 * _tmp285 +
                         _tmp241 * _tmp286 * _tmp99 - _tmp241 * _tmp306 + _tmp250 * _tmp288 +
                         _tmp257 * _tmp304 + _tmp264 * _tmp287 - _tmp308 * _tmp309 -
                         _tmp309 * _tmp310;
  const Scalar _tmp312 = std::pow(_tmp289, Scalar(-2));
  const Scalar _tmp313 =
      std::pow(Scalar(std::pow(_tmp296, Scalar(2)) * _tmp312 + 1), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp314 = _tmp296 * _tmp312;
  const Scalar _tmp315 =
      _tmp313 *
      (_tmp290 * (_tmp211 * _tmp294 + _tmp219 * _tmp291 + _tmp234 * _tmp293 + _tmp246 * _tmp292) -
       _tmp311 * _tmp314);
  const Scalar _tmp316 = Scalar(1.0) * std::cosh(_tmp302);
  const Scalar _tmp317 = std::cosh(_tmp301);
  const Scalar _tmp318 = Scalar(0.71007031138673404) * _tmp299 * _tmp312;
  const Scalar _tmp319 = Scalar(1.4083112389913199) * _tmp289;
  const Scalar _tmp320 = _tmp119 * _tmp143;
  const Scalar _tmp321 = _tmp119 * _tmp129;
  const Scalar _tmp322 = _tmp125 * _tmp321 + _tmp141 * _tmp320 + _tmp286 + _tmp307;
  const Scalar _tmp323 = Scalar(1.0) / (_tmp322);
  const Scalar _tmp324 = _tmp129 * _tmp62;
  const Scalar _tmp325 = _tmp143 * _tmp62;
  const Scalar _tmp326 = -_tmp124 * _tmp324 - _tmp140 * _tmp325 - _tmp145 * _tmp148 +
                         _tmp152 * _tmp154 + _tmp295 - _tmp48 * _tmp66;
  const Scalar _tmp327 = std::asinh(_tmp323 * _tmp326);
  const Scalar _tmp328 = Scalar(1.0) * _tmp327;
  const Scalar _tmp329 = Scalar(1.4083112389913199) * _tmp322;
  const Scalar _tmp330 =
      -_tmp327 * _tmp329 - std::sqrt(Scalar(std::pow(Scalar(-_tmp80 + p_c(0, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp82 + p_c(1, 0)), Scalar(2))));
  const Scalar _tmp331 = Scalar(0.71007031138673404) * _tmp323;
  const Scalar _tmp332 = _tmp330 * _tmp331;
  const Scalar _tmp333 = -std::sinh(_tmp328) - std::sinh(_tmp332);
  const Scalar _tmp334 = _tmp125 * _tmp129;
  const Scalar _tmp335 = -_tmp143 * _tmp233 - _tmp149 * _tmp208 - _tmp217 * _tmp334 +
                         _tmp231 * _tmp320 + _tmp245 * _tmp321 - _tmp305;
  const Scalar _tmp336 = Scalar(1.4083112389913199) * _tmp335;
  const Scalar _tmp337 = std::cosh(_tmp332);
  const Scalar _tmp338 = std::pow(_tmp322, Scalar(-2));
  const Scalar _tmp339 = Scalar(0.71007031138673404) * _tmp330 * _tmp338;
  const Scalar _tmp340 =
      std::pow(Scalar(std::pow(_tmp326, Scalar(2)) * _tmp338 + 1), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp341 = _tmp326 * _tmp338;
  const Scalar _tmp342 =
      _tmp340 *
      (_tmp323 * (_tmp148 * _tmp213 - _tmp148 * _tmp214 - _tmp227 * _tmp325 + _tmp229 * _tmp325 -
                  _tmp230 * _tmp325 + _tmp236 * _tmp324 - _tmp237 * _tmp324 - _tmp244 * _tmp324) -
       _tmp335 * _tmp341);
  const Scalar _tmp343 = Scalar(1.0) * std::cosh(_tmp328);
  const Scalar _tmp344 = _tmp102 * _tmp68;
  const Scalar _tmp345 = _tmp344 * _tmp95;
  const Scalar _tmp346 = _tmp185 * _tmp68;
  const Scalar _tmp347 = _tmp252 * _tmp346;
  const Scalar _tmp348 = _tmp346 * _tmp95;
  const Scalar _tmp349 = _tmp348 * _tmp99;
  const Scalar _tmp350 = _tmp113 * _tmp344;
  const Scalar _tmp351 = _tmp110 * _tmp68;
  const Scalar _tmp352 =
      -_tmp107 * _tmp172 + _tmp107 * _tmp85 - _tmp109 * _tmp349 + _tmp351 * _tmp85;
  const Scalar _tmp353 = -_tmp113 * _tmp349 + _tmp350 * _tmp85 - _tmp352 * _tmp59;
  const Scalar _tmp354 = _tmp216 * _tmp353;
  const Scalar _tmp355 = _tmp344 * _tmp96;
  const Scalar _tmp356 = _tmp157 * (-_tmp201 * _tmp203 + _tmp202 * _tmp71 - _tmp75);
  const Scalar _tmp357 = -_tmp100 * _tmp195 + _tmp356 * _tmp95;
  const Scalar _tmp358 = _tmp348 * _tmp96;
  const Scalar _tmp359 = -_tmp103 * _tmp357 - _tmp195 * _tmp98 + _tmp355 * _tmp85 +
                         _tmp356 * _tmp85 - _tmp358 * _tmp99;
  const Scalar _tmp360 = _tmp207 * _tmp359;
  const Scalar _tmp361 = _tmp360 * _tmp99;
  const Scalar _tmp362 = _tmp122 * _tmp353;
  const Scalar _tmp363 = _tmp139 * _tmp362;
  const Scalar _tmp364 = _tmp104 * _tmp354;
  const Scalar _tmp365 = -_tmp107 + _tmp158 * _tmp189 + _tmp188 - _tmp351;
  const Scalar _tmp366 = _tmp192 - _tmp350 - _tmp365 * _tmp59;
  const Scalar _tmp367 =
      _tmp123 * (_tmp137 * _tmp364 - _tmp138 * _tmp366 + _tmp195 * _tmp97 + _tmp198 +
                 _tmp222 * _tmp357 - _tmp223 * _tmp359 - _tmp355 - _tmp356);
  const Scalar _tmp368 = _tmp228 * _tmp360;
  const Scalar _tmp369 = _tmp363 + _tmp366 + _tmp367 - _tmp368;
  const Scalar _tmp370 = -_tmp159 * _tmp369 + _tmp168 * _tmp68 + _tmp263 * _tmp353 + _tmp68;
  const Scalar _tmp371 = _tmp232 * _tmp353;
  const Scalar _tmp372 = _tmp107 * _tmp238;
  const Scalar _tmp373 = _tmp240 * _tmp348;
  const Scalar _tmp374 = -Scalar(1.0) * _tmp113 * _tmp348 + _tmp372 * _tmp59 + _tmp373 * _tmp59;
  const Scalar _tmp375 = _tmp123 * (-_tmp114 * _tmp357 + _tmp115 * _tmp364 - _tmp120 * _tmp359 -
                                    _tmp138 * _tmp374 - Scalar(1.0) * _tmp358);
  const Scalar _tmp376 = _tmp121 * _tmp362;
  const Scalar _tmp377 = _tmp235 * _tmp360;
  const Scalar _tmp378 = _tmp374 + _tmp375 + _tmp376 - _tmp377;
  const Scalar _tmp379 = -_tmp159 * _tmp378 + _tmp166 * _tmp68 + _tmp249 * _tmp354;
  const Scalar _tmp380 = _tmp164 * _tmp346;
  const Scalar _tmp381 =
      _tmp248 * _tmp97 -
      _tmp255 * (_tmp158 * _tmp379 + _tmp165 * _tmp378 - _tmp251 * _tmp354 + _tmp252 * _tmp380) +
      _tmp256 * _tmp97 -
      _tmp261 * (_tmp171 * _tmp345 + _tmp238 * _tmp361 - _tmp258 * _tmp360 - _tmp259 * _tmp347) +
      _tmp262 * _tmp97 -
      _tmp266 * (_tmp158 * _tmp370 + _tmp165 * _tmp369 + _tmp265 * _tmp346 - _tmp371 * _tmp85) +
      _tmp267 * _tmp97 -
      _tmp270 * (-_tmp160 * _tmp347 + _tmp162 * _tmp345 + _tmp268 * _tmp354 - _tmp269 * _tmp354);
  const Scalar _tmp382 = _tmp119 * _tmp352;
  const Scalar _tmp383 = _tmp117 * _tmp353;
  const Scalar _tmp384 = -_tmp126 * _tmp369 + _tmp127 * _tmp363 + _tmp127 * _tmp367 -
                         _tmp127 * _tmp368 - _tmp141 * _tmp382 + _tmp232 * _tmp383 + _tmp365;
  const Scalar _tmp385 = _tmp216 * _tmp383;
  const Scalar _tmp386 = -_tmp125 * _tmp382 + _tmp125 * _tmp385 - _tmp126 * _tmp378 +
                         _tmp127 * _tmp375 + _tmp127 * _tmp376 - _tmp127 * _tmp377 - _tmp372 -
                         _tmp373;
  const Scalar _tmp387 =
      _tmp133 * _tmp362 - _tmp146 * _tmp352 - _tmp209 * _tmp360 + _tmp210 * _tmp360;
  const Scalar _tmp388 = _tmp212 * _tmp360;
  const Scalar _tmp389 = _tmp362 * _tmp63;
  const Scalar _tmp390 = -_tmp151 * _tmp382 + _tmp151 * _tmp385;
  const Scalar _tmp391 =
      _tmp180 *
      (_tmp175 *
           (_tmp130 * (_tmp375 * _tmp63 + _tmp376 * _tmp63 - _tmp377 * _tmp63 - _tmp386 * _tmp65) +
            _tmp144 * (_tmp363 * _tmp63 + _tmp367 * _tmp63 - _tmp368 * _tmp63 - _tmp384 * _tmp65) +
            _tmp149 * (-_tmp387 * _tmp65 - _tmp388 + _tmp389) - _tmp220 * _tmp390) -
       _tmp272 * _tmp381);
  const Scalar _tmp392 = Scalar(1.4083112389913199) * _tmp381;
  const Scalar _tmp393 = _tmp284 * _tmp354;
  const Scalar _tmp394 = _tmp103 * _tmp393 + _tmp129 * _tmp380 * _tmp95 - _tmp285 * _tmp348 -
                         _tmp286 * _tmp349 + _tmp287 * _tmp370 + _tmp288 * _tmp379 +
                         _tmp304 * _tmp361 + _tmp306 * _tmp348 + _tmp308 * _tmp344 +
                         _tmp310 * _tmp344;
  const Scalar _tmp395 = Scalar(1.4083112389913199) * _tmp394;
  const Scalar _tmp396 =
      _tmp313 *
      (_tmp290 * (_tmp291 * _tmp390 + _tmp292 * _tmp386 + _tmp293 * _tmp384 + _tmp294 * _tmp387) -
       _tmp314 * _tmp394);
  const Scalar _tmp397 = -_tmp143 * _tmp371 - _tmp149 * _tmp360 + _tmp320 * _tmp369 +
                         _tmp321 * _tmp378 - _tmp334 * _tmp354 - _tmp393;
  const Scalar _tmp398 = Scalar(1.4083112389913199) * _tmp397;
  const Scalar _tmp399 =
      _tmp340 *
      (_tmp323 * (_tmp148 * _tmp388 - _tmp148 * _tmp389 - _tmp324 * _tmp375 - _tmp324 * _tmp376 +
                  _tmp324 * _tmp377 - _tmp325 * _tmp363 - _tmp325 * _tmp367 + _tmp325 * _tmp368) -
       _tmp341 * _tmp397);

  // Output terms (1)
  Eigen::Matrix<Scalar, 4, 3> _res;

  _res(0, 0) = 0;
  _res(1, 0) =
      -_tmp18 * _tmp37 +
      _tmp275 * (-_tmp178 * _tmp273 -
                 _tmp282 * (-_tmp271 * _tmp277 + _tmp280 * (-_tmp176 * _tmp278 - _tmp18 * _tmp279 -
                                                            _tmp273 * _tmp275))) +
      _tmp278 * _tmp283;
  _res(2, 0) =
      Scalar(1.4083112389913199) * _tmp303 * _tmp311 +
      _tmp319 *
          (-_tmp315 * _tmp316 -
           _tmp317 * (_tmp300 * (-_tmp298 * _tmp311 - _tmp315 * _tmp319) - _tmp311 * _tmp318));
  _res(3, 0) =
      _tmp329 *
          (-_tmp337 * (_tmp331 * (-_tmp327 * _tmp336 - _tmp329 * _tmp342) - _tmp335 * _tmp339) -
           _tmp342 * _tmp343) +
      _tmp333 * _tmp336;
  _res(0, 1) = 0;
  _res(1, 1) =
      _tmp275 * (-_tmp178 * _tmp391 -
                 _tmp282 * (-_tmp277 * _tmp381 + _tmp280 * (-_tmp176 * _tmp392 - _tmp275 * _tmp391 -
                                                            _tmp279 * _tmp35))) +
      _tmp283 * _tmp392 - _tmp35 * _tmp37;
  _res(2, 1) =
      _tmp303 * _tmp395 + _tmp319 * (-_tmp316 * _tmp396 -
                                     _tmp317 * (_tmp300 * (-_tmp297 * _tmp395 - _tmp319 * _tmp396) -
                                                _tmp318 * _tmp394));
  _res(3, 1) =
      _tmp329 *
          (-_tmp337 * (_tmp331 * (-_tmp327 * _tmp398 - _tmp329 * _tmp399) - _tmp339 * _tmp397) -
           _tmp343 * _tmp399) +
      _tmp333 * _tmp398;
  _res(0, 2) = 0;
  _res(1, 2) = -_tmp27 * _tmp37;
  _res(2, 2) = 0;
  _res(3, 2) = 0;

  return _res;
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym

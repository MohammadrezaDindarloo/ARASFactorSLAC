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
 * Symbolic function: IK_residual_func_cost1_wrt_pc_Nl8
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
Eigen::Matrix<Scalar, 4, 3> IkResidualFuncCost1WrtPcNl8(
    const Scalar fh1, const Scalar fv1, const sym::Rot3<Scalar>& DeltaRot,
    const Eigen::Matrix<Scalar, 3, 1>& position_vector, const Eigen::Matrix<Scalar, 4, 1>& encoder,
    const Eigen::Matrix<Scalar, 3, 1>& p_a, const Eigen::Matrix<Scalar, 3, 1>& p_b,
    const Eigen::Matrix<Scalar, 3, 1>& p_c, const Eigen::Matrix<Scalar, 3, 1>& p_d,
    const sym::Rot3<Scalar>& Rot_init, const Scalar epsilon) {
  // Total ops: 1231

  // Unused inputs
  (void)encoder;
  (void)epsilon;

  // Input arrays
  const Eigen::Matrix<Scalar, 4, 1>& _DeltaRot = DeltaRot.Data();
  const Eigen::Matrix<Scalar, 4, 1>& _Rot_init = Rot_init.Data();

  // Intermediate terms (400)
  const Scalar _tmp0 = _DeltaRot[0] * _Rot_init[2] + _DeltaRot[1] * _Rot_init[3] -
                       _DeltaRot[2] * _Rot_init[0] + _DeltaRot[3] * _Rot_init[1];
  const Scalar _tmp1 = -2 * std::pow(_tmp0, Scalar(2));
  const Scalar _tmp2 = -_DeltaRot[0] * _Rot_init[1] + _DeltaRot[1] * _Rot_init[0] +
                       _DeltaRot[2] * _Rot_init[3] + _DeltaRot[3] * _Rot_init[2];
  const Scalar _tmp3 = 1 - 2 * std::pow(_tmp2, Scalar(2));
  const Scalar _tmp4 = Scalar(0.20999999999999999) * _tmp1 + Scalar(0.20999999999999999) * _tmp3;
  const Scalar _tmp5 = _DeltaRot[0] * _Rot_init[3] - _DeltaRot[1] * _Rot_init[2] +
                       _DeltaRot[2] * _Rot_init[1] + _DeltaRot[3] * _Rot_init[0];
  const Scalar _tmp6 = 2 * _tmp2 * _tmp5;
  const Scalar _tmp7 = -2 * _DeltaRot[0] * _Rot_init[0] - 2 * _DeltaRot[1] * _Rot_init[1] -
                       2 * _DeltaRot[2] * _Rot_init[2] + 2 * _DeltaRot[3] * _Rot_init[3];
  const Scalar _tmp8 = _tmp0 * _tmp7;
  const Scalar _tmp9 = _tmp6 + _tmp8;
  const Scalar _tmp10 = -Scalar(0.010999999999999999) * _tmp9;
  const Scalar _tmp11 = 2 * _tmp0;
  const Scalar _tmp12 = _tmp11 * _tmp5;
  const Scalar _tmp13 = _tmp2 * _tmp7;
  const Scalar _tmp14 = Scalar(0.20999999999999999) * _tmp12 - Scalar(0.20999999999999999) * _tmp13;
  const Scalar _tmp15 = _tmp10 + _tmp14;
  const Scalar _tmp16 = _tmp15 + _tmp4;
  const Scalar _tmp17 = _tmp16 + position_vector(0, 0);
  const Scalar _tmp18 = _tmp17 - p_c(0, 0);
  const Scalar _tmp19 = Scalar(1.0) / (_tmp18);
  const Scalar _tmp20 = std::pow(_tmp18, Scalar(2));
  const Scalar _tmp21 = -2 * std::pow(_tmp5, Scalar(2));
  const Scalar _tmp22 = Scalar(0.20999999999999999) * _tmp21 + Scalar(0.20999999999999999) * _tmp3;
  const Scalar _tmp23 = _tmp11 * _tmp2;
  const Scalar _tmp24 = _tmp5 * _tmp7;
  const Scalar _tmp25 = _tmp23 - _tmp24;
  const Scalar _tmp26 = -Scalar(0.010999999999999999) * _tmp25;
  const Scalar _tmp27 = Scalar(0.20999999999999999) * _tmp12 + Scalar(0.20999999999999999) * _tmp13;
  const Scalar _tmp28 = _tmp26 + _tmp27;
  const Scalar _tmp29 = _tmp22 + _tmp28;
  const Scalar _tmp30 = _tmp29 + position_vector(1, 0);
  const Scalar _tmp31 = _tmp30 - p_c(1, 0);
  const Scalar _tmp32 = std::pow(_tmp31, Scalar(2));
  const Scalar _tmp33 = _tmp20 + _tmp32;
  const Scalar _tmp34 = std::sqrt(_tmp33);
  const Scalar _tmp35 = _tmp19 * _tmp34;
  const Scalar _tmp36 = Scalar(1.0) / (_tmp34);
  const Scalar _tmp37 = _tmp29 * _tmp36;
  const Scalar _tmp38 = _tmp16 * _tmp36;
  const Scalar _tmp39 = -_tmp18 * _tmp37 + _tmp31 * _tmp38;
  const Scalar _tmp40 = _tmp35 * _tmp39;
  const Scalar _tmp41 = -_tmp4;
  const Scalar _tmp42 = _tmp15 + _tmp41;
  const Scalar _tmp43 = _tmp42 + position_vector(0, 0);
  const Scalar _tmp44 = _tmp43 - p_d(0, 0);
  const Scalar _tmp45 = _tmp26 - _tmp27;
  const Scalar _tmp46 = _tmp22 + _tmp45;
  const Scalar _tmp47 = _tmp46 + position_vector(1, 0);
  const Scalar _tmp48 = _tmp47 - p_d(1, 0);
  const Scalar _tmp49 = std::pow(Scalar(std::pow(_tmp44, Scalar(2)) + std::pow(_tmp48, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp50 = _tmp44 * _tmp49;
  const Scalar _tmp51 = _tmp48 * _tmp49;
  const Scalar _tmp52 = _tmp10 - _tmp14;
  const Scalar _tmp53 = _tmp41 + _tmp52;
  const Scalar _tmp54 = _tmp53 + position_vector(0, 0);
  const Scalar _tmp55 = _tmp54 - p_a(0, 0);
  const Scalar _tmp56 = std::pow(_tmp55, Scalar(2));
  const Scalar _tmp57 = -_tmp22;
  const Scalar _tmp58 = _tmp45 + _tmp57;
  const Scalar _tmp59 = _tmp58 + position_vector(1, 0);
  const Scalar _tmp60 = _tmp59 - p_a(1, 0);
  const Scalar _tmp61 = _tmp56 + std::pow(_tmp60, Scalar(2));
  const Scalar _tmp62 = std::pow(_tmp61, Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp63 = _tmp55 * _tmp62;
  const Scalar _tmp64 = _tmp60 * _tmp62;
  const Scalar _tmp65 = _tmp40 * _tmp63 - _tmp53 * _tmp64 + _tmp58 * _tmp63;
  const Scalar _tmp66 = _tmp19 * _tmp31;
  const Scalar _tmp67 = _tmp63 * _tmp66;
  const Scalar _tmp68 = -_tmp64 + _tmp67;
  const Scalar _tmp69 = Scalar(1.0) / (_tmp68);
  const Scalar _tmp70 = _tmp50 * _tmp66 - _tmp51;
  const Scalar _tmp71 = _tmp69 * _tmp70;
  const Scalar _tmp72 = _tmp40 * _tmp50 - _tmp42 * _tmp51 + _tmp46 * _tmp50 - _tmp65 * _tmp71;
  const Scalar _tmp73 = Scalar(1.0) / (_tmp72);
  const Scalar _tmp74 = Scalar(1.0) * _tmp73;
  const Scalar _tmp75 = Scalar(1.0) * _tmp69;
  const Scalar _tmp76 = _tmp70 * _tmp75;
  const Scalar _tmp77 = _tmp73 * _tmp76;
  const Scalar _tmp78 = _tmp4 + _tmp52;
  const Scalar _tmp79 = _tmp78 - p_b(0, 0) + position_vector(0, 0);
  const Scalar _tmp80 = _tmp28 + _tmp57;
  const Scalar _tmp81 = _tmp80 - p_b(1, 0) + position_vector(1, 0);
  const Scalar _tmp82 = std::pow(Scalar(std::pow(_tmp79, Scalar(2)) + std::pow(_tmp81, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp83 = _tmp81 * _tmp82;
  const Scalar _tmp84 = _tmp79 * _tmp82;
  const Scalar _tmp85 = fh1 * (_tmp78 * _tmp83 - _tmp80 * _tmp84);
  const Scalar _tmp86 = _tmp85 * (_tmp50 * _tmp74 - _tmp63 * _tmp77);
  const Scalar _tmp87 = Scalar(0.20999999999999999) * _tmp23 + Scalar(0.20999999999999999) * _tmp24;
  const Scalar _tmp88 = -Scalar(0.010999999999999999) * _tmp1 -
                        Scalar(0.010999999999999999) * _tmp21 + Scalar(-0.010999999999999999);
  const Scalar _tmp89 = Scalar(0.20999999999999999) * _tmp6 - Scalar(0.20999999999999999) * _tmp8;
  const Scalar _tmp90 = _tmp88 + _tmp89;
  const Scalar _tmp91 = _tmp87 + _tmp90;
  const Scalar _tmp92 = _tmp50 * _tmp91;
  const Scalar _tmp93 = _tmp88 - _tmp89;
  const Scalar _tmp94 = _tmp87 + _tmp93;
  const Scalar _tmp95 = _tmp66 * _tmp91;
  const Scalar _tmp96 = -_tmp87;
  const Scalar _tmp97 = _tmp93 + _tmp96;
  const Scalar _tmp98 = -_tmp63 * _tmp95 + _tmp64 * _tmp97;
  const Scalar _tmp99 = _tmp69 * _tmp98;
  const Scalar _tmp100 = _tmp51 * _tmp94 - _tmp66 * _tmp92 - _tmp70 * _tmp99;
  const Scalar _tmp101 = Scalar(1.0) * _tmp29;
  const Scalar _tmp102 = -_tmp101;
  const Scalar _tmp103 = Scalar(1.0) / (_tmp102 + _tmp58);
  const Scalar _tmp104 = Scalar(1.0) * _tmp16;
  const Scalar _tmp105 = _tmp104 - _tmp53;
  const Scalar _tmp106 = _tmp103 * _tmp105;
  const Scalar _tmp107 = _tmp63 * _tmp91;
  const Scalar _tmp108 = _tmp107 - _tmp63 * _tmp97;
  const Scalar _tmp109 = -_tmp100 * _tmp106 - _tmp108 * _tmp71 - _tmp50 * _tmp94 + _tmp92;
  const Scalar _tmp110 = Scalar(1.0) / (_tmp109);
  const Scalar _tmp111 = Scalar(1.0) * _tmp103;
  const Scalar _tmp112 = _tmp105 * _tmp111;
  const Scalar _tmp113 = -_tmp108 * _tmp75 + _tmp112 * _tmp99;
  const Scalar _tmp114 = _tmp110 * _tmp113;
  const Scalar _tmp115 = -_tmp114 * _tmp72 - _tmp65 * _tmp75;
  const Scalar _tmp116 = _tmp109 * _tmp73;
  const Scalar _tmp117 = _tmp115 * _tmp116;
  const Scalar _tmp118 = _tmp113 + _tmp117;
  const Scalar _tmp119 = _tmp110 * _tmp70;
  const Scalar _tmp120 = -_tmp118 * _tmp119 + Scalar(1.0);
  const Scalar _tmp121 = _tmp63 * _tmp69;
  const Scalar _tmp122 = _tmp110 * _tmp118;
  const Scalar _tmp123 = _tmp122 * _tmp50;
  const Scalar _tmp124 = _tmp83 * fh1;
  const Scalar _tmp125 = _tmp124 * (_tmp120 * _tmp121 + _tmp123);
  const Scalar _tmp126 = _tmp101 * _tmp106 + _tmp104;
  const Scalar _tmp127 = 0;
  const Scalar _tmp128 = _tmp110 * _tmp127;
  const Scalar _tmp129 = _tmp63 * _tmp71;
  const Scalar _tmp130 = _tmp128 * _tmp50;
  const Scalar _tmp131 = Scalar(333.54000000000002) - fv1;
  const Scalar _tmp132 = _tmp131 * (-_tmp128 * _tmp129 + _tmp130);
  const Scalar _tmp133 = _tmp66 * _tmp69;
  const Scalar _tmp134 = _tmp66 * _tmp99 + _tmp95;
  const Scalar _tmp135 = -_tmp106 * _tmp134 + _tmp108 * _tmp133 - _tmp91;
  const Scalar _tmp136 = _tmp110 * _tmp135;
  const Scalar _tmp137 = _tmp65 * _tmp69;
  const Scalar _tmp138 = -_tmp136 * _tmp72 + _tmp137 * _tmp66 - _tmp40;
  const Scalar _tmp139 = _tmp116 * _tmp138;
  const Scalar _tmp140 = _tmp135 + _tmp139;
  const Scalar _tmp141 = _tmp110 * _tmp140;
  const Scalar _tmp142 = -_tmp141 * _tmp70 - _tmp66;
  const Scalar _tmp143 = _tmp141 * _tmp50;
  const Scalar _tmp144 = _tmp84 * fh1;
  const Scalar _tmp145 = _tmp144 * (_tmp121 * _tmp142 + _tmp143 + Scalar(1.0));
  const Scalar _tmp146 = -_tmp125 * _tmp35 - _tmp132 * _tmp35 - _tmp145 * _tmp35 - _tmp35 * _tmp86;
  const Scalar _tmp147 = Scalar(1.0) / (_tmp146);
  const Scalar _tmp148 = Scalar(0.71007031138673404) * _tmp147;
  const Scalar _tmp149 = _tmp100 * _tmp110;
  const Scalar _tmp150 = _tmp102 + _tmp46;
  const Scalar _tmp151 = _tmp106 * _tmp150;
  const Scalar _tmp152 = Scalar(1.0) / (_tmp104 - _tmp151 - _tmp42);
  const Scalar _tmp153 = _tmp126 * _tmp152;
  const Scalar _tmp154 = _tmp102 - _tmp127 * _tmp149 - _tmp150 * _tmp153;
  const Scalar _tmp155 = Scalar(1.0) * _tmp152;
  const Scalar _tmp156 = _tmp150 * _tmp152;
  const Scalar _tmp157 = _tmp117 * _tmp156 - _tmp118 * _tmp149 - Scalar(1.0) * _tmp99;
  const Scalar _tmp158 = Scalar(1.0) * _tmp124;
  const Scalar _tmp159 = _tmp116 * _tmp155;
  const Scalar _tmp160 = _tmp150 * _tmp155;
  const Scalar _tmp161 = -_tmp100 * _tmp74 + _tmp116 * _tmp160;
  const Scalar _tmp162 = Scalar(1.0) * _tmp85;
  const Scalar _tmp163 = fh1 * (_tmp90 + _tmp96);
  const Scalar _tmp164 = _tmp163 * _tmp84 + _tmp78 * fv1 + Scalar(40.024799999999999) * _tmp9;
  const Scalar _tmp165 = _tmp134 + _tmp139 * _tmp156 - _tmp140 * _tmp149;
  const Scalar _tmp166 = Scalar(1.0) * _tmp144;
  const Scalar _tmp167 = -_tmp163 * _tmp83 - Scalar(40.024799999999999) * _tmp25 - _tmp80 * fv1;
  const Scalar _tmp168 = _tmp151 * _tmp155 + Scalar(1.0);
  const Scalar _tmp169 = _tmp106 * _tmp155;
  const Scalar _tmp170 =
      Scalar(1.0) * _tmp131 * (-_tmp111 * _tmp154 - _tmp126 * _tmp155 + Scalar(1.0)) +
      _tmp158 * (-_tmp111 * _tmp157 + _tmp117 * _tmp155) +
      _tmp162 * (-_tmp111 * _tmp161 + _tmp159) +
      Scalar(1.0) * _tmp164 * (_tmp103 * _tmp160 - _tmp155) +
      _tmp166 * (-_tmp111 * _tmp165 + _tmp139 * _tmp155) +
      Scalar(1.0) * _tmp167 * (-_tmp111 * _tmp168 + _tmp169);
  const Scalar _tmp171 = std::asinh(_tmp147 * _tmp170);
  const Scalar _tmp172 = Scalar(1.4083112389913199) * _tmp146;
  const Scalar _tmp173 = -_tmp30 + p_c(1, 0);
  const Scalar _tmp174 = -_tmp17 + p_c(0, 0);
  const Scalar _tmp175 =
      std::sqrt(Scalar(std::pow(_tmp173, Scalar(2)) + std::pow(_tmp174, Scalar(2))));
  const Scalar _tmp176 = -_tmp171 * _tmp172 - _tmp175;
  const Scalar _tmp177 = _tmp148 * _tmp176;
  const Scalar _tmp178 = Scalar(1.0) * _tmp171;
  const Scalar _tmp179 = _tmp148 * p_c(2, 0) - std::cosh(_tmp177) + std::cosh(_tmp178);
  const Scalar _tmp180 = Scalar(1.0) / (_tmp20);
  const Scalar _tmp181 = _tmp180 * _tmp34;
  const Scalar _tmp182 = _tmp132 * _tmp36;
  const Scalar _tmp183 = _tmp180 * _tmp31;
  const Scalar _tmp184 = _tmp108 * _tmp183;
  const Scalar _tmp185 = _tmp184 * _tmp69;
  const Scalar _tmp186 = _tmp183 * _tmp99;
  const Scalar _tmp187 = _tmp183 * _tmp91;
  const Scalar _tmp188 = std::pow(_tmp68, Scalar(-2));
  const Scalar _tmp189 = _tmp188 * _tmp63;
  const Scalar _tmp190 = _tmp183 * _tmp189;
  const Scalar _tmp191 = _tmp190 * _tmp98;
  const Scalar _tmp192 = _tmp129 * _tmp187 - _tmp183 * _tmp92 - _tmp186 * _tmp50 + _tmp191 * _tmp70;
  const Scalar _tmp193 = _tmp184 * _tmp189;
  const Scalar _tmp194 = -_tmp106 * _tmp192 - _tmp185 * _tmp50 + _tmp193 * _tmp70;
  const Scalar _tmp195 = _tmp115 * _tmp73;
  const Scalar _tmp196 = _tmp194 * _tmp195;
  const Scalar _tmp197 = _tmp121 * _tmp187;
  const Scalar _tmp198 = -_tmp112 * _tmp191 - _tmp112 * _tmp197 + Scalar(1.0) * _tmp193;
  const Scalar _tmp199 = _tmp110 * _tmp72;
  const Scalar _tmp200 = std::pow(_tmp109, Scalar(-2));
  const Scalar _tmp201 = _tmp200 * _tmp72;
  const Scalar _tmp202 = _tmp113 * _tmp201;
  const Scalar _tmp203 = _tmp181 * _tmp39;
  const Scalar _tmp204 = _tmp137 * _tmp183;
  const Scalar _tmp205 = _tmp36 * _tmp39;
  const Scalar _tmp206 = _tmp205 * _tmp50;
  const Scalar _tmp207 = _tmp190 * _tmp65;
  const Scalar _tmp208 = std::pow(_tmp33, Scalar(Scalar(-3) / Scalar(2)));
  const Scalar _tmp209 = _tmp208 * _tmp29;
  const Scalar _tmp210 = _tmp16 * _tmp208;
  const Scalar _tmp211 = _tmp18 * _tmp31;
  const Scalar _tmp212 = _tmp35 * (-_tmp20 * _tmp209 + _tmp210 * _tmp211 + _tmp37);
  const Scalar _tmp213 = _tmp203 * _tmp63 - _tmp205 * _tmp63 + _tmp212 * _tmp63;
  const Scalar _tmp214 = _tmp203 * _tmp50 - _tmp204 * _tmp50 - _tmp206 + _tmp207 * _tmp70 +
                         _tmp212 * _tmp50 - _tmp213 * _tmp71;
  const Scalar _tmp215 = _tmp116 * (-_tmp114 * _tmp214 + _tmp194 * _tmp202 - _tmp198 * _tmp199 +
                                    Scalar(1.0) * _tmp207 - _tmp213 * _tmp75);
  const Scalar _tmp216 = std::pow(_tmp72, Scalar(-2));
  const Scalar _tmp217 = _tmp109 * _tmp216;
  const Scalar _tmp218 = _tmp115 * _tmp217;
  const Scalar _tmp219 = _tmp214 * _tmp218;
  const Scalar _tmp220 = _tmp196 + _tmp198 + _tmp215 - _tmp219;
  const Scalar _tmp221 = _tmp110 * _tmp50;
  const Scalar _tmp222 = _tmp118 * _tmp200;
  const Scalar _tmp223 = _tmp222 * _tmp70;
  const Scalar _tmp224 = -_tmp119 * _tmp220 - _tmp123 * _tmp183 + _tmp194 * _tmp223;
  const Scalar _tmp225 = _tmp194 * _tmp50;
  const Scalar _tmp226 = _tmp56 / _tmp61;
  const Scalar _tmp227 = _tmp120 * _tmp226;
  const Scalar _tmp228 = _tmp183 * _tmp188;
  const Scalar _tmp229 = _tmp124 * _tmp35;
  const Scalar _tmp230 = _tmp140 * _tmp200;
  const Scalar _tmp231 = _tmp226 * _tmp228;
  const Scalar _tmp232 = _tmp138 * _tmp73;
  const Scalar _tmp233 = _tmp194 * _tmp232;
  const Scalar _tmp234 = _tmp135 * _tmp201;
  const Scalar _tmp235 = _tmp32 / [&]() {
    const Scalar base = _tmp18;
    return base * base * base;
  }();
  const Scalar _tmp236 = _tmp188 * _tmp235;
  const Scalar _tmp237 = _tmp63 * _tmp98;
  const Scalar _tmp238 = -_tmp107 * _tmp235 * _tmp69 + _tmp186 + _tmp187 - _tmp236 * _tmp237;
  const Scalar _tmp239 = _tmp108 * _tmp63;
  const Scalar _tmp240 = -_tmp106 * _tmp238 + _tmp185 - _tmp236 * _tmp239;
  const Scalar _tmp241 = _tmp63 * _tmp65;
  const Scalar _tmp242 =
      _tmp116 * (_tmp133 * _tmp213 - _tmp136 * _tmp214 + _tmp194 * _tmp234 - _tmp199 * _tmp240 -
                 _tmp203 + _tmp204 + _tmp205 - _tmp212 - _tmp236 * _tmp241);
  const Scalar _tmp243 = _tmp138 * _tmp217;
  const Scalar _tmp244 = _tmp214 * _tmp243;
  const Scalar _tmp245 = _tmp233 + _tmp240 + _tmp242 - _tmp244;
  const Scalar _tmp246 = _tmp230 * _tmp70;
  const Scalar _tmp247 = -_tmp119 * _tmp245 - _tmp143 * _tmp183 - _tmp183 + _tmp194 * _tmp246;
  const Scalar _tmp248 = _tmp144 * _tmp35;
  const Scalar _tmp249 = _tmp63 * _tmp75;
  const Scalar _tmp250 = _tmp50 * _tmp73;
  const Scalar _tmp251 = _tmp183 * _tmp250;
  const Scalar _tmp252 = _tmp216 * _tmp76;
  const Scalar _tmp253 = _tmp214 * _tmp252;
  const Scalar _tmp254 = Scalar(1.0) * _tmp216;
  const Scalar _tmp255 = _tmp254 * _tmp50;
  const Scalar _tmp256 = _tmp35 * _tmp85;
  const Scalar _tmp257 = _tmp127 * _tmp200;
  const Scalar _tmp258 = _tmp129 * _tmp257;
  const Scalar _tmp259 = _tmp119 * _tmp127;
  const Scalar _tmp260 = _tmp131 * _tmp35;
  const Scalar _tmp261 = _tmp36 * _tmp86;
  const Scalar _tmp262 = _tmp125 * _tmp36;
  const Scalar _tmp263 = _tmp145 * _tmp36;
  const Scalar _tmp264 =
      -_tmp125 * _tmp181 - _tmp132 * _tmp181 - _tmp145 * _tmp181 - _tmp181 * _tmp86 + _tmp182 -
      _tmp229 * (_tmp121 * _tmp224 + _tmp220 * _tmp221 - _tmp222 * _tmp225 - _tmp227 * _tmp228) -
      _tmp248 * (_tmp121 * _tmp247 - _tmp142 * _tmp231 + _tmp221 * _tmp245 - _tmp225 * _tmp230) -
      _tmp256 *
          (-_tmp214 * _tmp255 + _tmp231 * _tmp70 * _tmp74 - _tmp249 * _tmp251 + _tmp253 * _tmp63) -
      _tmp260 * (-_tmp121 * _tmp130 * _tmp183 + _tmp194 * _tmp258 - _tmp225 * _tmp257 +
                 _tmp231 * _tmp259) +
      _tmp261 + _tmp262 + _tmp263;
  const Scalar _tmp265 = Scalar(1.4083112389913199) * _tmp264;
  const Scalar _tmp266 = std::pow(_tmp146, Scalar(-2));
  const Scalar _tmp267 = _tmp170 * _tmp266;
  const Scalar _tmp268 = _tmp100 * _tmp194;
  const Scalar _tmp269 = -_tmp141 * _tmp192 - _tmp149 * _tmp245 + _tmp156 * _tmp233 +
                         _tmp156 * _tmp242 - _tmp156 * _tmp244 + _tmp230 * _tmp268 + _tmp238;
  const Scalar _tmp270 = _tmp155 * _tmp73;
  const Scalar _tmp271 = _tmp194 * _tmp270;
  const Scalar _tmp272 = _tmp155 * _tmp217;
  const Scalar _tmp273 = _tmp214 * _tmp272;
  const Scalar _tmp274 = _tmp160 * _tmp217;
  const Scalar _tmp275 = _tmp100 * _tmp254;
  const Scalar _tmp276 =
      _tmp150 * _tmp271 - _tmp192 * _tmp74 - _tmp214 * _tmp274 + _tmp214 * _tmp275;
  const Scalar _tmp277 = -_tmp128 * _tmp192 + _tmp257 * _tmp268;
  const Scalar _tmp278 = _tmp111 * _tmp131;
  const Scalar _tmp279 = -_tmp122 * _tmp192 - _tmp149 * _tmp220 + _tmp156 * _tmp196 +
                         _tmp156 * _tmp215 - _tmp156 * _tmp219 + _tmp187 * _tmp249 +
                         Scalar(1.0) * _tmp191 + _tmp222 * _tmp268;
  const Scalar _tmp280 = _tmp147 * (_tmp158 * (-_tmp111 * _tmp279 + _tmp155 * _tmp196 +
                                               _tmp155 * _tmp215 - _tmp155 * _tmp219) +
                                    _tmp162 * (-_tmp111 * _tmp276 + _tmp271 - _tmp273) +
                                    _tmp166 * (-_tmp111 * _tmp269 + _tmp155 * _tmp233 +
                                               _tmp155 * _tmp242 - _tmp155 * _tmp244) -
                                    _tmp277 * _tmp278) -
                         _tmp264 * _tmp267;
  const Scalar _tmp281 =
      std::pow(Scalar(std::pow(_tmp170, Scalar(2)) * _tmp266 + 1), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp282 = Scalar(1.0) * _tmp281 * std::sinh(_tmp178);
  const Scalar _tmp283 = Scalar(0.71007031138673404) * _tmp266;
  const Scalar _tmp284 = _tmp283 * p_c(2, 0);
  const Scalar _tmp285 = _tmp176 * _tmp283;
  const Scalar _tmp286 = _tmp172 * _tmp281;
  const Scalar _tmp287 = Scalar(1.0) / (_tmp175);
  const Scalar _tmp288 = std::sinh(_tmp177);
  const Scalar _tmp289 = _tmp144 * _tmp69;
  const Scalar _tmp290 = _tmp124 * _tmp69;
  const Scalar _tmp291 = _tmp128 * _tmp131;
  const Scalar _tmp292 = _tmp120 * _tmp290 + _tmp142 * _tmp289 - _tmp291 * _tmp71 - _tmp77 * _tmp85;
  const Scalar _tmp293 = Scalar(1.0) / (_tmp292);
  const Scalar _tmp294 = Scalar(0.71007031138673404) * _tmp293;
  const Scalar _tmp295 = _tmp103 * _tmp85;
  const Scalar _tmp296 = _tmp155 * _tmp164;
  const Scalar _tmp297 = _tmp103 * _tmp131;
  const Scalar _tmp298 = _tmp103 * _tmp124;
  const Scalar _tmp299 = _tmp103 * _tmp144;
  const Scalar _tmp300 = -_tmp103 * _tmp150 * _tmp296 + _tmp103 * _tmp167 * _tmp168 +
                         _tmp154 * _tmp297 + _tmp157 * _tmp298 + _tmp161 * _tmp295 +
                         _tmp165 * _tmp299;
  const Scalar _tmp301 = std::asinh(_tmp293 * _tmp300);
  const Scalar _tmp302 = Scalar(1.4083112389913199) * _tmp292;
  const Scalar _tmp303 =
      -_tmp301 * _tmp302 - std::sqrt(Scalar(std::pow(Scalar(-_tmp54 + p_a(0, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp59 + p_a(1, 0)), Scalar(2))));
  const Scalar _tmp304 = _tmp294 * _tmp303;
  const Scalar _tmp305 = Scalar(1.0) * _tmp301;
  const Scalar _tmp306 = _tmp294 * p_a(2, 0) - std::cosh(_tmp304) + std::cosh(_tmp305);
  const Scalar _tmp307 = _tmp291 * _tmp50;
  const Scalar _tmp308 = _tmp131 * _tmp257;
  const Scalar _tmp309 = _tmp194 * _tmp308;
  const Scalar _tmp310 = _tmp120 * _tmp124;
  const Scalar _tmp311 = _tmp74 * _tmp85;
  const Scalar _tmp312 = _tmp75 * _tmp85;
  const Scalar _tmp313 = _tmp131 * _tmp190 * _tmp259 - _tmp142 * _tmp144 * _tmp190 -
                         _tmp183 * _tmp307 * _tmp69 - _tmp190 * _tmp310 +
                         _tmp190 * _tmp311 * _tmp70 + _tmp224 * _tmp290 + _tmp247 * _tmp289 -
                         _tmp251 * _tmp312 + _tmp253 * _tmp85 + _tmp309 * _tmp71;
  const Scalar _tmp314 = Scalar(1.4083112389913199) * _tmp313;
  const Scalar _tmp315 = Scalar(1.0) * std::sinh(_tmp305);
  const Scalar _tmp316 = std::pow(_tmp292, Scalar(-2));
  const Scalar _tmp317 =
      std::pow(Scalar(std::pow(_tmp300, Scalar(2)) * _tmp316 + 1), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp318 = _tmp300 * _tmp316;
  const Scalar _tmp319 =
      _tmp317 *
      (_tmp293 * (_tmp269 * _tmp299 + _tmp276 * _tmp295 + _tmp277 * _tmp297 + _tmp279 * _tmp298) -
       _tmp313 * _tmp318);
  const Scalar _tmp320 = Scalar(0.71007031138673404) * _tmp316;
  const Scalar _tmp321 = _tmp320 * p_a(2, 0);
  const Scalar _tmp322 = _tmp303 * _tmp320;
  const Scalar _tmp323 = std::sinh(_tmp304);
  const Scalar _tmp324 = _tmp124 * _tmp222;
  const Scalar _tmp325 = _tmp110 * _tmp124;
  const Scalar _tmp326 = _tmp162 * _tmp216;
  const Scalar _tmp327 = _tmp110 * _tmp144;
  const Scalar _tmp328 = _tmp144 * _tmp230;
  const Scalar _tmp329 = -_tmp194 * _tmp324 - _tmp194 * _tmp328 - _tmp214 * _tmp326 +
                         _tmp220 * _tmp325 + _tmp245 * _tmp327 - _tmp309;
  const Scalar _tmp330 = _tmp124 * _tmp152;
  const Scalar _tmp331 = _tmp144 * _tmp152;
  const Scalar _tmp332 = -_tmp117 * _tmp330 + _tmp131 * _tmp153 - _tmp139 * _tmp331 -
                         _tmp159 * _tmp85 - _tmp167 * _tmp169 + _tmp296;
  const Scalar _tmp333 = _tmp122 * _tmp124 + _tmp141 * _tmp144 + _tmp291 + _tmp311;
  const Scalar _tmp334 = Scalar(1.0) / (_tmp333);
  const Scalar _tmp335 = std::asinh(_tmp332 * _tmp334);
  const Scalar _tmp336 = Scalar(1.4083112389913199) * _tmp333;
  const Scalar _tmp337 =
      -_tmp335 * _tmp336 - std::sqrt(Scalar(std::pow(Scalar(-_tmp43 + p_d(0, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp47 + p_d(1, 0)), Scalar(2))));
  const Scalar _tmp338 = Scalar(0.71007031138673404) * _tmp334;
  const Scalar _tmp339 = _tmp337 * _tmp338;
  const Scalar _tmp340 = Scalar(1.0) * _tmp335;
  const Scalar _tmp341 = Scalar(1.4083112389913199) * _tmp338 * p_d(2, 0) -
                         Scalar(1.4083112389913199) * std::cosh(_tmp339) +
                         Scalar(1.4083112389913199) * std::cosh(_tmp340);
  const Scalar _tmp342 = Scalar(1.0) * std::sinh(_tmp340);
  const Scalar _tmp343 = std::pow(_tmp333, Scalar(-2));
  const Scalar _tmp344 =
      std::pow(Scalar(std::pow(_tmp332, Scalar(2)) * _tmp343 + 1), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp345 = _tmp332 * _tmp343;
  const Scalar _tmp346 =
      _tmp344 *
      (-_tmp329 * _tmp345 +
       _tmp334 * (-_tmp196 * _tmp330 - _tmp215 * _tmp330 + _tmp219 * _tmp330 - _tmp233 * _tmp331 -
                  _tmp242 * _tmp331 + _tmp244 * _tmp331 - _tmp271 * _tmp85 + _tmp273 * _tmp85));
  const Scalar _tmp347 = std::sinh(_tmp339);
  const Scalar _tmp348 = Scalar(1.4083112389913199) * _tmp335;
  const Scalar _tmp349 = Scalar(0.71007031138673404) * _tmp343;
  const Scalar _tmp350 = _tmp337 * _tmp349;
  const Scalar _tmp351 = _tmp349 * p_d(2, 0);
  const Scalar _tmp352 = _tmp19 * _tmp99;
  const Scalar _tmp353 = _tmp19 * _tmp91;
  const Scalar _tmp354 = _tmp191 + _tmp197 - _tmp352 - _tmp353;
  const Scalar _tmp355 = _tmp19 * _tmp69;
  const Scalar _tmp356 = _tmp108 * _tmp355;
  const Scalar _tmp357 = -_tmp106 * _tmp354 + _tmp193 - _tmp356;
  const Scalar _tmp358 = _tmp35 * (-_tmp209 * _tmp211 + _tmp210 * _tmp32 - _tmp38);
  const Scalar _tmp359 = -_tmp205 * _tmp67 + _tmp358 * _tmp63;
  const Scalar _tmp360 = _tmp355 * _tmp65;
  const Scalar _tmp361 = _tmp188 * _tmp19;
  const Scalar _tmp362 = _tmp361 * _tmp70;
  const Scalar _tmp363 = -_tmp206 * _tmp66 - _tmp241 * _tmp362 + _tmp358 * _tmp50 -
                         _tmp359 * _tmp71 + _tmp360 * _tmp50;
  const Scalar _tmp364 = _tmp239 * _tmp361;
  const Scalar _tmp365 =
      -_tmp129 * _tmp353 - _tmp237 * _tmp362 + _tmp352 * _tmp50 + _tmp353 * _tmp50;
  const Scalar _tmp366 = -_tmp106 * _tmp365 + _tmp356 * _tmp50 - _tmp364 * _tmp70;
  const Scalar _tmp367 =
      _tmp116 * (_tmp133 * _tmp359 - _tmp136 * _tmp363 - _tmp199 * _tmp357 + _tmp205 * _tmp66 +
                 _tmp207 + _tmp234 * _tmp366 - _tmp358 - _tmp360);
  const Scalar _tmp368 = _tmp243 * _tmp363;
  const Scalar _tmp369 = _tmp232 * _tmp366;
  const Scalar _tmp370 = _tmp357 + _tmp367 - _tmp368 + _tmp369;
  const Scalar _tmp371 = -_tmp119 * _tmp370 + _tmp143 * _tmp19 + _tmp19 + _tmp246 * _tmp366;
  const Scalar _tmp372 = _tmp366 * _tmp50;
  const Scalar _tmp373 = _tmp142 * _tmp361;
  const Scalar _tmp374 = _tmp259 * _tmp361;
  const Scalar _tmp375 = _tmp218 * _tmp363;
  const Scalar _tmp376 = _tmp195 * _tmp366;
  const Scalar _tmp377 =
      _tmp112 * _tmp121 * _tmp353 + _tmp112 * _tmp237 * _tmp361 - Scalar(1.0) * _tmp364;
  const Scalar _tmp378 = Scalar(1.0) * _tmp361;
  const Scalar _tmp379 = _tmp116 * (-_tmp114 * _tmp363 - _tmp199 * _tmp377 + _tmp202 * _tmp366 -
                                    _tmp241 * _tmp378 - _tmp359 * _tmp75);
  const Scalar _tmp380 = -_tmp375 + _tmp376 + _tmp377 + _tmp379;
  const Scalar _tmp381 = -_tmp119 * _tmp380 + _tmp123 * _tmp19 + _tmp223 * _tmp366;
  const Scalar _tmp382 = _tmp19 * _tmp250;
  const Scalar _tmp383 = _tmp252 * _tmp363;
  const Scalar _tmp384 =
      _tmp182 * _tmp66 -
      _tmp229 * (_tmp121 * _tmp381 + _tmp221 * _tmp380 - _tmp222 * _tmp372 + _tmp227 * _tmp361) -
      _tmp248 * (_tmp121 * _tmp371 + _tmp221 * _tmp370 + _tmp226 * _tmp373 - _tmp230 * _tmp372) -
      _tmp256 *
          (-_tmp226 * _tmp362 * _tmp74 + _tmp249 * _tmp382 - _tmp255 * _tmp363 + _tmp383 * _tmp63) -
      _tmp260 *
          (_tmp130 * _tmp355 * _tmp63 - _tmp226 * _tmp374 - _tmp257 * _tmp372 + _tmp258 * _tmp366) +
      _tmp261 * _tmp66 + _tmp262 * _tmp66 + _tmp263 * _tmp66;
  const Scalar _tmp385 = _tmp100 * _tmp366;
  const Scalar _tmp386 = -_tmp128 * _tmp365 + _tmp257 * _tmp385;
  const Scalar _tmp387 = _tmp270 * _tmp366;
  const Scalar _tmp388 =
      _tmp150 * _tmp387 - _tmp274 * _tmp363 + _tmp275 * _tmp363 - _tmp365 * _tmp74;
  const Scalar _tmp389 = _tmp272 * _tmp363;
  const Scalar _tmp390 = -_tmp141 * _tmp365 - _tmp149 * _tmp370 + _tmp156 * _tmp367 -
                         _tmp156 * _tmp368 + _tmp156 * _tmp369 + _tmp230 * _tmp385 + _tmp354;
  const Scalar _tmp391 = -_tmp122 * _tmp365 - _tmp149 * _tmp380 - _tmp156 * _tmp375 +
                         _tmp156 * _tmp376 + _tmp156 * _tmp379 + _tmp222 * _tmp385 -
                         _tmp237 * _tmp378 - _tmp249 * _tmp353;
  const Scalar _tmp392 = _tmp147 * (_tmp158 * (-_tmp111 * _tmp391 - _tmp155 * _tmp375 +
                                               _tmp155 * _tmp376 + _tmp155 * _tmp379) +
                                    _tmp162 * (-_tmp111 * _tmp388 + _tmp387 - _tmp389) +
                                    _tmp166 * (-_tmp111 * _tmp390 + _tmp155 * _tmp367 -
                                               _tmp155 * _tmp368 + _tmp155 * _tmp369) -
                                    _tmp278 * _tmp386) -
                         _tmp267 * _tmp384;
  const Scalar _tmp393 = Scalar(1.4083112389913199) * _tmp384;
  const Scalar _tmp394 = _tmp308 * _tmp366;
  const Scalar _tmp395 = -_tmp131 * _tmp374 * _tmp63 + _tmp144 * _tmp373 * _tmp63 +
                         _tmp289 * _tmp371 + _tmp290 * _tmp381 + _tmp307 * _tmp355 +
                         _tmp310 * _tmp361 * _tmp63 - _tmp311 * _tmp362 * _tmp63 +
                         _tmp312 * _tmp382 + _tmp383 * _tmp85 + _tmp394 * _tmp71;
  const Scalar _tmp396 = Scalar(1.4083112389913199) * _tmp395;
  const Scalar _tmp397 =
      _tmp317 *
      (_tmp293 * (_tmp295 * _tmp388 + _tmp297 * _tmp386 + _tmp298 * _tmp391 + _tmp299 * _tmp390) -
       _tmp318 * _tmp395);
  const Scalar _tmp398 = -_tmp324 * _tmp366 + _tmp325 * _tmp380 - _tmp326 * _tmp363 +
                         _tmp327 * _tmp370 - _tmp328 * _tmp366 - _tmp394;
  const Scalar _tmp399 =
      _tmp344 *
      (_tmp334 * (_tmp330 * _tmp375 - _tmp330 * _tmp376 - _tmp330 * _tmp379 - _tmp331 * _tmp367 +
                  _tmp331 * _tmp368 - _tmp331 * _tmp369 - _tmp387 * _tmp85 + _tmp389 * _tmp85) -
       _tmp345 * _tmp398);

  // Output terms (1)
  Eigen::Matrix<Scalar, 4, 3> _res;

  _res(0, 0) = 0;
  _res(1, 0) =
      -_tmp172 *
          (-_tmp264 * _tmp284 + _tmp280 * _tmp282 -
           _tmp288 * (_tmp148 * (-_tmp171 * _tmp265 - _tmp174 * _tmp287 - _tmp280 * _tmp286) -
                      _tmp264 * _tmp285)) -
      _tmp179 * _tmp265;
  _res(2, 0) =
      -_tmp302 *
          (-_tmp313 * _tmp321 + _tmp315 * _tmp319 -
           _tmp323 * (_tmp294 * (-_tmp301 * _tmp314 - _tmp302 * _tmp319) - _tmp313 * _tmp322)) -
      _tmp306 * _tmp314;
  _res(3, 0) =
      -_tmp329 * _tmp341 -
      _tmp336 *
          (-_tmp329 * _tmp351 + _tmp342 * _tmp346 -
           _tmp347 * (-_tmp329 * _tmp350 + _tmp338 * (-_tmp329 * _tmp348 - _tmp336 * _tmp346)));
  _res(0, 1) = 0;
  _res(1, 1) =
      -_tmp172 *
          (_tmp282 * _tmp392 - _tmp284 * _tmp384 -
           _tmp288 * (_tmp148 * (-_tmp171 * _tmp393 - _tmp173 * _tmp287 - _tmp286 * _tmp392) -
                      _tmp285 * _tmp384)) -
      _tmp179 * _tmp393;
  _res(2, 1) =
      -_tmp302 *
          (_tmp315 * _tmp397 - _tmp321 * _tmp395 -
           _tmp323 * (_tmp294 * (-_tmp301 * _tmp396 - _tmp302 * _tmp397) - _tmp322 * _tmp395)) -
      _tmp306 * _tmp396;
  _res(3, 1) =
      -_tmp336 *
          (_tmp342 * _tmp399 -
           _tmp347 * (_tmp338 * (-_tmp336 * _tmp399 - _tmp348 * _tmp398) - _tmp350 * _tmp398) -
           _tmp351 * _tmp398) -
      _tmp341 * _tmp398;
  _res(0, 2) = 0;
  _res(1, 2) = Scalar(-1.0);
  _res(2, 2) = 0;
  _res(3, 2) = 0;

  return _res;
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym

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
 * Symbolic function: IK_residual_func_cost1_wrt_pc_Nl1
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
Eigen::Matrix<Scalar, 4, 3> IkResidualFuncCost1WrtPcNl1(
    const Scalar fh1, const Scalar fv1, const sym::Rot3<Scalar>& DeltaRot,
    const Eigen::Matrix<Scalar, 3, 1>& position_vector, const Eigen::Matrix<Scalar, 4, 1>& encoder,
    const Eigen::Matrix<Scalar, 3, 1>& p_a, const Eigen::Matrix<Scalar, 3, 1>& p_b,
    const Eigen::Matrix<Scalar, 3, 1>& p_c, const Eigen::Matrix<Scalar, 3, 1>& p_d,
    const sym::Rot3<Scalar>& Rot_init, const Scalar epsilon) {
  // Total ops: 1066

  // Unused inputs
  (void)encoder;
  (void)epsilon;

  // Input arrays
  const Eigen::Matrix<Scalar, 4, 1>& _DeltaRot = DeltaRot.Data();
  const Eigen::Matrix<Scalar, 4, 1>& _Rot_init = Rot_init.Data();

  // Intermediate terms (339)
  const Scalar _tmp0 = _DeltaRot[0] * _Rot_init[2] + _DeltaRot[1] * _Rot_init[3] -
                       _DeltaRot[2] * _Rot_init[0] + _DeltaRot[3] * _Rot_init[1];
  const Scalar _tmp1 = _DeltaRot[0] * _Rot_init[3] - _DeltaRot[1] * _Rot_init[2] +
                       _DeltaRot[2] * _Rot_init[1] + _DeltaRot[3] * _Rot_init[0];
  const Scalar _tmp2 = 2 * _tmp1;
  const Scalar _tmp3 = _tmp0 * _tmp2;
  const Scalar _tmp4 = -_DeltaRot[0] * _Rot_init[1] + _DeltaRot[1] * _Rot_init[0] +
                       _DeltaRot[2] * _Rot_init[3] + _DeltaRot[3] * _Rot_init[2];
  const Scalar _tmp5 = -_DeltaRot[0] * _Rot_init[0] - _DeltaRot[1] * _Rot_init[1] -
                       _DeltaRot[2] * _Rot_init[2] + _DeltaRot[3] * _Rot_init[3];
  const Scalar _tmp6 = 2 * _tmp5;
  const Scalar _tmp7 = _tmp4 * _tmp6;
  const Scalar _tmp8 = Scalar(0.20999999999999999) * _tmp3 - Scalar(0.20999999999999999) * _tmp7;
  const Scalar _tmp9 = -2 * std::pow(_tmp4, Scalar(2));
  const Scalar _tmp10 = 1 - 2 * std::pow(_tmp0, Scalar(2));
  const Scalar _tmp11 = Scalar(0.20999999999999999) * _tmp10 + Scalar(0.20999999999999999) * _tmp9;
  const Scalar _tmp12 = _tmp2 * _tmp4;
  const Scalar _tmp13 = _tmp0 * _tmp6;
  const Scalar _tmp14 = _tmp12 + _tmp13;
  const Scalar _tmp15 = -Scalar(0.010999999999999999) * _tmp14;
  const Scalar _tmp16 = _tmp11 + _tmp15;
  const Scalar _tmp17 = _tmp16 + _tmp8;
  const Scalar _tmp18 = _tmp17 + position_vector(0, 0);
  const Scalar _tmp19 = _tmp18 - p_c(0, 0);
  const Scalar _tmp20 = -_tmp8;
  const Scalar _tmp21 = _tmp16 + _tmp20;
  const Scalar _tmp22 = _tmp21 + position_vector(0, 0);
  const Scalar _tmp23 = _tmp22 - p_b(0, 0);
  const Scalar _tmp24 = Scalar(0.20999999999999999) * _tmp3 + Scalar(0.20999999999999999) * _tmp7;
  const Scalar _tmp25 = 2 * _tmp0 * _tmp4;
  const Scalar _tmp26 = _tmp2 * _tmp5;
  const Scalar _tmp27 = _tmp25 - _tmp26;
  const Scalar _tmp28 = -Scalar(0.010999999999999999) * _tmp27;
  const Scalar _tmp29 = -2 * std::pow(_tmp1, Scalar(2));
  const Scalar _tmp30 = Scalar(0.20999999999999999) * _tmp29 + Scalar(0.20999999999999999) * _tmp9 +
                        Scalar(0.20999999999999999);
  const Scalar _tmp31 = _tmp28 - _tmp30;
  const Scalar _tmp32 = _tmp24 + _tmp31;
  const Scalar _tmp33 = _tmp32 + position_vector(1, 0);
  const Scalar _tmp34 = _tmp33 - p_b(1, 0);
  const Scalar _tmp35 =
      std::sqrt(Scalar(std::pow(_tmp23, Scalar(2)) + std::pow(_tmp34, Scalar(2))));
  const Scalar _tmp36 = Scalar(1.0) / (_tmp35);
  const Scalar _tmp37 = Scalar(1.0) / (_tmp23);
  const Scalar _tmp38 = _tmp35 * _tmp37;
  const Scalar _tmp39 = _tmp38 * (_tmp21 * _tmp34 * _tmp36 - _tmp23 * _tmp32 * _tmp36);
  const Scalar _tmp40 = Scalar(0.20999999999999999) * _tmp25 + Scalar(0.20999999999999999) * _tmp26;
  const Scalar _tmp41 = -_tmp40;
  const Scalar _tmp42 =
      -Scalar(0.010999999999999999) * _tmp10 - Scalar(0.010999999999999999) * _tmp29;
  const Scalar _tmp43 = Scalar(0.20999999999999999) * _tmp12 - Scalar(0.20999999999999999) * _tmp13;
  const Scalar _tmp44 = _tmp42 + _tmp43;
  const Scalar _tmp45 = _tmp41 + _tmp44;
  const Scalar _tmp46 = _tmp34 * _tmp37;
  const Scalar _tmp47 = _tmp45 * _tmp46;
  const Scalar _tmp48 = -_tmp11 + _tmp15;
  const Scalar _tmp49 = _tmp48 + _tmp8;
  const Scalar _tmp50 = _tmp49 + position_vector(0, 0);
  const Scalar _tmp51 = _tmp50 - p_d(0, 0);
  const Scalar _tmp52 = -_tmp24;
  const Scalar _tmp53 = _tmp28 + _tmp30;
  const Scalar _tmp54 = _tmp52 + _tmp53;
  const Scalar _tmp55 = _tmp54 + position_vector(1, 0);
  const Scalar _tmp56 = _tmp55 - p_d(1, 0);
  const Scalar _tmp57 = std::pow(Scalar(std::pow(_tmp51, Scalar(2)) + std::pow(_tmp56, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp58 = _tmp51 * _tmp57;
  const Scalar _tmp59 = _tmp56 * _tmp57;
  const Scalar _tmp60 = Scalar(1.0) / (_tmp46 * _tmp58 - _tmp59);
  const Scalar _tmp61 = _tmp42 - _tmp43;
  const Scalar _tmp62 = _tmp40 + _tmp61;
  const Scalar _tmp63 = _tmp45 * _tmp58;
  const Scalar _tmp64 = -_tmp46 * _tmp63 + _tmp59 * _tmp62;
  const Scalar _tmp65 = _tmp60 * _tmp64;
  const Scalar _tmp66 = _tmp46 * _tmp65 + _tmp47;
  const Scalar _tmp67 = Scalar(1.0) * _tmp32;
  const Scalar _tmp68 = -_tmp67;
  const Scalar _tmp69 = Scalar(1.0) / (_tmp54 + _tmp68);
  const Scalar _tmp70 = Scalar(1.0) * _tmp21;
  const Scalar _tmp71 = -_tmp49 + _tmp70;
  const Scalar _tmp72 = _tmp69 * _tmp71;
  const Scalar _tmp73 = -_tmp58 * _tmp62 + _tmp63;
  const Scalar _tmp74 = _tmp60 * _tmp73;
  const Scalar _tmp75 = -_tmp45 + _tmp46 * _tmp74 - _tmp66 * _tmp72;
  const Scalar _tmp76 = std::pow(_tmp19, Scalar(2));
  const Scalar _tmp77 = _tmp24 + _tmp53;
  const Scalar _tmp78 = _tmp77 + position_vector(1, 0);
  const Scalar _tmp79 = _tmp78 - p_c(1, 0);
  const Scalar _tmp80 = std::pow(_tmp79, Scalar(2));
  const Scalar _tmp81 = _tmp76 + _tmp80;
  const Scalar _tmp82 = std::pow(_tmp81, Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp83 = _tmp40 + _tmp44;
  const Scalar _tmp84 = _tmp82 * _tmp83;
  const Scalar _tmp85 = _tmp46 * _tmp82;
  const Scalar _tmp86 = _tmp79 * _tmp82;
  const Scalar _tmp87 = _tmp19 * _tmp85 - _tmp86;
  const Scalar _tmp88 = _tmp45 * _tmp82;
  const Scalar _tmp89 = _tmp19 * _tmp88;
  const Scalar _tmp90 = -_tmp46 * _tmp89 - _tmp65 * _tmp87 + _tmp83 * _tmp86;
  const Scalar _tmp91 = -_tmp19 * _tmp84 - _tmp72 * _tmp90 - _tmp74 * _tmp87 + _tmp89;
  const Scalar _tmp92 = Scalar(1.0) / (_tmp91);
  const Scalar _tmp93 = _tmp77 * _tmp82;
  const Scalar _tmp94 = _tmp17 * _tmp82;
  const Scalar _tmp95 = _tmp39 * _tmp82;
  const Scalar _tmp96 = _tmp39 * _tmp58 - _tmp49 * _tmp59 + _tmp54 * _tmp58;
  const Scalar _tmp97 = _tmp60 * _tmp96;
  const Scalar _tmp98 = _tmp19 * _tmp93 + _tmp19 * _tmp95 - _tmp79 * _tmp94 - _tmp87 * _tmp97;
  const Scalar _tmp99 = _tmp92 * _tmp98;
  const Scalar _tmp100 = -_tmp39 + _tmp46 * _tmp97 - _tmp75 * _tmp99;
  const Scalar _tmp101 = Scalar(1.0) / (_tmp98);
  const Scalar _tmp102 = _tmp101 * _tmp91;
  const Scalar _tmp103 = _tmp100 * _tmp102;
  const Scalar _tmp104 = _tmp103 + _tmp75;
  const Scalar _tmp105 = _tmp82 * _tmp92;
  const Scalar _tmp106 = _tmp104 * _tmp105;
  const Scalar _tmp107 = _tmp87 * _tmp92;
  const Scalar _tmp108 = -_tmp104 * _tmp107 - _tmp46;
  const Scalar _tmp109 = _tmp58 * _tmp60;
  const Scalar _tmp110 = _tmp20 + _tmp48;
  const Scalar _tmp111 = _tmp110 - p_a(0, 0) + position_vector(0, 0);
  const Scalar _tmp112 = _tmp31 + _tmp52;
  const Scalar _tmp113 = _tmp112 - p_a(1, 0) + position_vector(1, 0);
  const Scalar _tmp114 =
      std::pow(Scalar(std::pow(_tmp111, Scalar(2)) + std::pow(_tmp113, Scalar(2))),
               Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp115 = _tmp111 * _tmp114;
  const Scalar _tmp116 = _tmp115 * fh1;
  const Scalar _tmp117 = _tmp116 * _tmp38;
  const Scalar _tmp118 = Scalar(1.0) * _tmp101;
  const Scalar _tmp119 = _tmp118 * _tmp82;
  const Scalar _tmp120 = Scalar(1.0) * _tmp60;
  const Scalar _tmp121 = _tmp101 * _tmp120;
  const Scalar _tmp122 = _tmp121 * _tmp87;
  const Scalar _tmp123 = _tmp113 * _tmp114;
  const Scalar _tmp124 = fh1 * (_tmp110 * _tmp123 - _tmp112 * _tmp115);
  const Scalar _tmp125 = _tmp124 * _tmp38;
  const Scalar _tmp126 = _tmp120 * _tmp64;
  const Scalar _tmp127 = -_tmp120 * _tmp73 + _tmp126 * _tmp72;
  const Scalar _tmp128 = -_tmp120 * _tmp96 - _tmp127 * _tmp99;
  const Scalar _tmp129 = _tmp102 * _tmp128;
  const Scalar _tmp130 = _tmp127 + _tmp129;
  const Scalar _tmp131 = _tmp130 * _tmp92;
  const Scalar _tmp132 = -_tmp131 * _tmp87 + Scalar(1.0);
  const Scalar _tmp133 = _tmp105 * _tmp130;
  const Scalar _tmp134 = _tmp123 * fh1;
  const Scalar _tmp135 = _tmp134 * _tmp38;
  const Scalar _tmp136 = _tmp67 * _tmp72 + _tmp70;
  const Scalar _tmp137 = 0;
  const Scalar _tmp138 = _tmp109 * _tmp137;
  const Scalar _tmp139 = _tmp105 * _tmp137;
  const Scalar _tmp140 = Scalar(333.54000000000002) - fv1;
  const Scalar _tmp141 = _tmp140 * _tmp38;
  const Scalar _tmp142 = -_tmp117 * (_tmp106 * _tmp19 + _tmp108 * _tmp109 + Scalar(1.0)) -
                         _tmp125 * (_tmp119 * _tmp19 - _tmp122 * _tmp58) -
                         _tmp135 * (_tmp109 * _tmp132 + _tmp133 * _tmp19) -
                         _tmp141 * (-_tmp107 * _tmp138 + _tmp139 * _tmp19);
  const Scalar _tmp143 = std::pow(_tmp142, Scalar(-2));
  const Scalar _tmp144 = _tmp68 + _tmp77;
  const Scalar _tmp145 = _tmp144 * _tmp72;
  const Scalar _tmp146 = Scalar(1.0) / (-_tmp145 - _tmp17 + _tmp70);
  const Scalar _tmp147 = Scalar(1.0) * _tmp146;
  const Scalar _tmp148 = _tmp102 * _tmp147;
  const Scalar _tmp149 = -_tmp118 * _tmp90 + _tmp144 * _tmp148;
  const Scalar _tmp150 = Scalar(1.0) * _tmp69;
  const Scalar _tmp151 = Scalar(1.0) * _tmp124;
  const Scalar _tmp152 = fh1 * (_tmp41 + _tmp61);
  const Scalar _tmp153 = -_tmp112 * fv1 - _tmp123 * _tmp152 - Scalar(40.024799999999999) * _tmp27;
  const Scalar _tmp154 = _tmp145 * _tmp147 + Scalar(1.0);
  const Scalar _tmp155 = _tmp110 * fv1 + _tmp115 * _tmp152 + Scalar(40.024799999999999) * _tmp14;
  const Scalar _tmp156 = _tmp144 * _tmp69;
  const Scalar _tmp157 = _tmp90 * _tmp92;
  const Scalar _tmp158 = _tmp144 * _tmp146;
  const Scalar _tmp159 = -_tmp126 + _tmp129 * _tmp158 - _tmp130 * _tmp157;
  const Scalar _tmp160 = Scalar(1.0) * fh1;
  const Scalar _tmp161 = _tmp123 * _tmp160;
  const Scalar _tmp162 = _tmp103 * _tmp158 - _tmp104 * _tmp157 + _tmp66;
  const Scalar _tmp163 = _tmp115 * _tmp160;
  const Scalar _tmp164 = _tmp136 * _tmp146;
  const Scalar _tmp165 = -_tmp137 * _tmp157 - _tmp144 * _tmp164 + _tmp68;
  const Scalar _tmp166 =
      Scalar(1.0) * _tmp140 * (-_tmp136 * _tmp147 - _tmp150 * _tmp165 + Scalar(1.0)) +
      _tmp151 * (_tmp148 - _tmp149 * _tmp150) +
      Scalar(1.0) * _tmp153 * (_tmp147 * _tmp72 - _tmp150 * _tmp154) +
      Scalar(1.0) * _tmp155 * (_tmp147 * _tmp156 - _tmp147) +
      _tmp161 * (_tmp129 * _tmp147 - _tmp150 * _tmp159) +
      _tmp163 * (_tmp103 * _tmp147 - _tmp150 * _tmp162);
  const Scalar _tmp167 =
      std::pow(Scalar(_tmp143 * std::pow(_tmp166, Scalar(2)) + 1), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp168 = Scalar(1.0) / (_tmp142);
  const Scalar _tmp169 = std::pow(_tmp81, Scalar(Scalar(-3) / Scalar(2)));
  const Scalar _tmp170 = _tmp169 * _tmp76;
  const Scalar _tmp171 = _tmp169 * _tmp19 * _tmp79;
  const Scalar _tmp172 = _tmp170 * _tmp46 - _tmp171 - _tmp85;
  const Scalar _tmp173 = _tmp171 * _tmp83;
  const Scalar _tmp174 = _tmp170 * _tmp45;
  const Scalar _tmp175 = -_tmp172 * _tmp65 + _tmp173 - _tmp174 * _tmp46 + _tmp46 * _tmp88;
  const Scalar _tmp176 =
      -_tmp170 * _tmp83 - _tmp172 * _tmp74 + _tmp174 - _tmp175 * _tmp72 + _tmp84 - _tmp88;
  const Scalar _tmp177 = _tmp101 * _tmp176;
  const Scalar _tmp178 = _tmp128 * _tmp177;
  const Scalar _tmp179 =
      -_tmp17 * _tmp171 + _tmp170 * _tmp39 + _tmp170 * _tmp77 - _tmp172 * _tmp97 - _tmp93 - _tmp95;
  const Scalar _tmp180 = std::pow(_tmp98, Scalar(-2));
  const Scalar _tmp181 = _tmp180 * _tmp91;
  const Scalar _tmp182 = _tmp179 * _tmp181;
  const Scalar _tmp183 = _tmp128 * _tmp182;
  const Scalar _tmp184 = _tmp127 * _tmp92;
  const Scalar _tmp185 = _tmp127 * _tmp98;
  const Scalar _tmp186 = std::pow(_tmp91, Scalar(-2));
  const Scalar _tmp187 = _tmp176 * _tmp186;
  const Scalar _tmp188 = _tmp102 * (-_tmp179 * _tmp184 + _tmp185 * _tmp187);
  const Scalar _tmp189 = _tmp178 - _tmp183 + _tmp188;
  const Scalar _tmp190 = _tmp186 * _tmp90;
  const Scalar _tmp191 = _tmp176 * _tmp190;
  const Scalar _tmp192 = _tmp130 * _tmp191 - _tmp131 * _tmp175 - _tmp157 * _tmp189 +
                         _tmp158 * _tmp178 - _tmp158 * _tmp183 + _tmp158 * _tmp188;
  const Scalar _tmp193 = _tmp75 * _tmp92;
  const Scalar _tmp194 = _tmp75 * _tmp98;
  const Scalar _tmp195 = _tmp102 * (-_tmp179 * _tmp193 + _tmp187 * _tmp194);
  const Scalar _tmp196 = _tmp100 * _tmp177;
  const Scalar _tmp197 = _tmp100 * _tmp182;
  const Scalar _tmp198 = _tmp195 + _tmp196 - _tmp197;
  const Scalar _tmp199 = _tmp104 * _tmp176;
  const Scalar _tmp200 = _tmp104 * _tmp92;
  const Scalar _tmp201 = -_tmp157 * _tmp198 + _tmp158 * _tmp195 + _tmp158 * _tmp196 -
                         _tmp158 * _tmp197 - _tmp175 * _tmp200 + _tmp190 * _tmp199;
  const Scalar _tmp202 = _tmp147 * _tmp182;
  const Scalar _tmp203 = Scalar(1.0) * _tmp180;
  const Scalar _tmp204 = _tmp203 * _tmp90;
  const Scalar _tmp205 = _tmp147 * _tmp177;
  const Scalar _tmp206 =
      -_tmp118 * _tmp175 - _tmp144 * _tmp202 + _tmp144 * _tmp205 + _tmp179 * _tmp204;
  const Scalar _tmp207 = _tmp137 * _tmp92;
  const Scalar _tmp208 = _tmp140 * _tmp69;
  const Scalar _tmp209 = _tmp208 * (_tmp137 * _tmp191 - _tmp175 * _tmp207);
  const Scalar _tmp210 = _tmp19 * _tmp82;
  const Scalar _tmp211 = _tmp109 * _tmp207;
  const Scalar _tmp212 = _tmp186 * _tmp87;
  const Scalar _tmp213 = _tmp176 * _tmp212;
  const Scalar _tmp214 = _tmp105 * _tmp19;
  const Scalar _tmp215 = _tmp130 * _tmp187;
  const Scalar _tmp216 = -_tmp107 * _tmp189 + _tmp130 * _tmp213 - _tmp131 * _tmp172;
  const Scalar _tmp217 = _tmp121 * _tmp172;
  const Scalar _tmp218 = _tmp120 * _tmp180 * _tmp87;
  const Scalar _tmp219 = _tmp179 * _tmp218;
  const Scalar _tmp220 = _tmp203 * _tmp210;
  const Scalar _tmp221 = -_tmp107 * _tmp198 - _tmp172 * _tmp200 + _tmp199 * _tmp212;
  const Scalar _tmp222 = _tmp186 * _tmp199;
  const Scalar _tmp223 = -_tmp117 * (-_tmp106 + _tmp109 * _tmp221 + _tmp170 * _tmp200 +
                                     _tmp198 * _tmp214 - _tmp210 * _tmp222) -
                         _tmp125 * (_tmp118 * _tmp170 - _tmp119 - _tmp179 * _tmp220 -
                                    _tmp217 * _tmp58 + _tmp219 * _tmp58) -
                         _tmp135 * (_tmp109 * _tmp216 + _tmp131 * _tmp170 - _tmp133 +
                                    _tmp189 * _tmp214 - _tmp210 * _tmp215) -
                         _tmp141 * (-_tmp137 * _tmp187 * _tmp210 + _tmp138 * _tmp213 - _tmp139 +
                                    _tmp170 * _tmp207 - _tmp172 * _tmp211);
  const Scalar _tmp224 = _tmp143 * _tmp166;
  const Scalar _tmp225 = _tmp167 * (_tmp168 * (_tmp151 * (-_tmp150 * _tmp206 - _tmp202 + _tmp205) +
                                               _tmp161 * (_tmp147 * _tmp178 - _tmp147 * _tmp183 +
                                                          _tmp147 * _tmp188 - _tmp150 * _tmp192) +
                                               _tmp163 * (_tmp147 * _tmp195 + _tmp147 * _tmp196 -
                                                          _tmp147 * _tmp197 - _tmp150 * _tmp201) -
                                               Scalar(1.0) * _tmp209) -
                                    _tmp223 * _tmp224);
  const Scalar _tmp226 = std::asinh(_tmp166 * _tmp168);
  const Scalar _tmp227 = Scalar(1.0) * _tmp226;
  const Scalar _tmp228 = Scalar(1.0) * std::sinh(_tmp227);
  const Scalar _tmp229 = Scalar(0.71007031138673404) * _tmp143;
  const Scalar _tmp230 = _tmp229 * p_b(2, 0);
  const Scalar _tmp231 = Scalar(1.4083112389913199) * _tmp226;
  const Scalar _tmp232 =
      -_tmp142 * _tmp231 - std::sqrt(Scalar(std::pow(Scalar(-_tmp22 + p_b(0, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp33 + p_b(1, 0)), Scalar(2))));
  const Scalar _tmp233 = Scalar(0.71007031138673404) * _tmp168;
  const Scalar _tmp234 = _tmp232 * _tmp233;
  const Scalar _tmp235 = std::sinh(_tmp234);
  const Scalar _tmp236 = _tmp229 * _tmp232;
  const Scalar _tmp237 = Scalar(1.4083112389913199) * _tmp142;
  const Scalar _tmp238 = _tmp233 * p_b(2, 0) + std::cosh(_tmp227) - std::cosh(_tmp234);
  const Scalar _tmp239 = _tmp140 * _tmp207;
  const Scalar _tmp240 = _tmp239 * _tmp60;
  const Scalar _tmp241 = _tmp137 * _tmp140;
  const Scalar _tmp242 = _tmp241 * _tmp60;
  const Scalar _tmp243 = _tmp60 * fh1;
  const Scalar _tmp244 = _tmp123 * _tmp243;
  const Scalar _tmp245 = _tmp115 * _tmp243;
  const Scalar _tmp246 = -_tmp124 * _tmp217 + _tmp124 * _tmp219 - _tmp172 * _tmp240 +
                         _tmp213 * _tmp242 + _tmp216 * _tmp244 + _tmp221 * _tmp245;
  const Scalar _tmp247 =
      -_tmp107 * _tmp242 + _tmp108 * _tmp245 - _tmp122 * _tmp124 + _tmp132 * _tmp244;
  const Scalar _tmp248 = Scalar(1.0) / (_tmp247);
  const Scalar _tmp249 = _tmp69 * fh1;
  const Scalar _tmp250 = _tmp123 * _tmp249;
  const Scalar _tmp251 = _tmp153 * _tmp69;
  const Scalar _tmp252 = _tmp115 * _tmp249;
  const Scalar _tmp253 = _tmp124 * _tmp69;
  const Scalar _tmp254 = _tmp147 * _tmp155;
  const Scalar _tmp255 = _tmp149 * _tmp253 + _tmp154 * _tmp251 - _tmp156 * _tmp254 +
                         _tmp159 * _tmp250 + _tmp162 * _tmp252 + _tmp165 * _tmp208;
  const Scalar _tmp256 = std::asinh(_tmp248 * _tmp255);
  const Scalar _tmp257 = Scalar(1.4083112389913199) * _tmp247;
  const Scalar _tmp258 =
      -_tmp256 * _tmp257 - std::sqrt(Scalar(std::pow(Scalar(-_tmp50 + p_d(0, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp55 + p_d(1, 0)), Scalar(2))));
  const Scalar _tmp259 = Scalar(0.71007031138673404) * _tmp248;
  const Scalar _tmp260 = _tmp258 * _tmp259;
  const Scalar _tmp261 = Scalar(1.0) * _tmp256;
  const Scalar _tmp262 = Scalar(1.4083112389913199) * _tmp259 * p_d(2, 0) -
                         Scalar(1.4083112389913199) * std::cosh(_tmp260) +
                         Scalar(1.4083112389913199) * std::cosh(_tmp261);
  const Scalar _tmp263 = std::sinh(_tmp260);
  const Scalar _tmp264 = std::pow(_tmp247, Scalar(-2));
  const Scalar _tmp265 = Scalar(0.71007031138673404) * _tmp264;
  const Scalar _tmp266 = _tmp258 * _tmp265;
  const Scalar _tmp267 = Scalar(1.4083112389913199) * _tmp256;
  const Scalar _tmp268 =
      std::pow(Scalar(std::pow(_tmp255, Scalar(2)) * _tmp264 + 1), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp269 = _tmp255 * _tmp264;
  const Scalar _tmp270 =
      _tmp268 * (-_tmp246 * _tmp269 +
                 _tmp248 * (_tmp192 * _tmp250 + _tmp201 * _tmp252 + _tmp206 * _tmp253 + _tmp209));
  const Scalar _tmp271 = Scalar(1.0) * std::sinh(_tmp261);
  const Scalar _tmp272 = _tmp265 * p_d(2, 0);
  const Scalar _tmp273 = _tmp151 * _tmp180;
  const Scalar _tmp274 = _tmp116 * _tmp92;
  const Scalar _tmp275 = _tmp134 * _tmp92;
  const Scalar _tmp276 = -_tmp116 * _tmp222 - _tmp134 * _tmp215 - _tmp179 * _tmp273 -
                         _tmp187 * _tmp241 + _tmp189 * _tmp275 + _tmp198 * _tmp274;
  const Scalar _tmp277 = _tmp116 * _tmp200 + _tmp118 * _tmp124 + _tmp131 * _tmp134 + _tmp239;
  const Scalar _tmp278 = Scalar(1.0) / (_tmp277);
  const Scalar _tmp279 = Scalar(0.71007031138673404) * _tmp278;
  const Scalar _tmp280 = _tmp116 * _tmp146;
  const Scalar _tmp281 = _tmp134 * _tmp146;
  const Scalar _tmp282 = -_tmp103 * _tmp280 - _tmp124 * _tmp148 - _tmp129 * _tmp281 +
                         _tmp140 * _tmp164 - _tmp147 * _tmp251 * _tmp71 + _tmp254;
  const Scalar _tmp283 = std::asinh(_tmp278 * _tmp282);
  const Scalar _tmp284 = Scalar(1.4083112389913199) * _tmp283;
  const Scalar _tmp285 = -_tmp78 + p_c(1, 0);
  const Scalar _tmp286 = -_tmp18 + p_c(0, 0);
  const Scalar _tmp287 =
      std::sqrt(Scalar(std::pow(_tmp285, Scalar(2)) + std::pow(_tmp286, Scalar(2))));
  const Scalar _tmp288 = -_tmp277 * _tmp284 - _tmp287;
  const Scalar _tmp289 = _tmp279 * _tmp288;
  const Scalar _tmp290 = Scalar(1.0) * _tmp283;
  const Scalar _tmp291 = Scalar(1.4083112389913199) * _tmp279 * p_c(2, 0) -
                         Scalar(1.4083112389913199) * std::cosh(_tmp289) +
                         Scalar(1.4083112389913199) * std::cosh(_tmp290);
  const Scalar _tmp292 = std::pow(_tmp277, Scalar(-2));
  const Scalar _tmp293 = Scalar(0.71007031138673404) * _tmp292;
  const Scalar _tmp294 = _tmp293 * p_c(2, 0);
  const Scalar _tmp295 = Scalar(1.4083112389913199) * _tmp277;
  const Scalar _tmp296 =
      std::pow(Scalar(std::pow(_tmp282, Scalar(2)) * _tmp292 + 1), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp297 = _tmp282 * _tmp292;
  const Scalar _tmp298 =
      _tmp296 *
      (-_tmp276 * _tmp297 +
       _tmp278 * (_tmp124 * _tmp202 - _tmp124 * _tmp205 - _tmp178 * _tmp281 + _tmp183 * _tmp281 -
                  _tmp188 * _tmp281 - _tmp195 * _tmp280 - _tmp196 * _tmp280 + _tmp197 * _tmp280));
  const Scalar _tmp299 = Scalar(1.0) / (_tmp287);
  const Scalar _tmp300 = _tmp288 * _tmp293;
  const Scalar _tmp301 = std::sinh(_tmp289);
  const Scalar _tmp302 = Scalar(1.0) * std::sinh(_tmp290);
  const Scalar _tmp303 = _tmp169 * _tmp80;
  const Scalar _tmp304 = _tmp171 * _tmp46 - _tmp303 + _tmp82;
  const Scalar _tmp305 =
      -_tmp17 * _tmp303 + _tmp171 * _tmp39 + _tmp171 * _tmp77 - _tmp304 * _tmp97 + _tmp94;
  const Scalar _tmp306 = -_tmp171 * _tmp47 + _tmp303 * _tmp83 - _tmp304 * _tmp65 - _tmp84;
  const Scalar _tmp307 = _tmp171 * _tmp45 - _tmp173 - _tmp304 * _tmp74 - _tmp306 * _tmp72;
  const Scalar _tmp308 = _tmp186 * _tmp307;
  const Scalar _tmp309 = _tmp102 * (-_tmp184 * _tmp305 + _tmp185 * _tmp308);
  const Scalar _tmp310 = _tmp101 * _tmp307;
  const Scalar _tmp311 = _tmp128 * _tmp310;
  const Scalar _tmp312 = _tmp181 * _tmp305;
  const Scalar _tmp313 = _tmp128 * _tmp312;
  const Scalar _tmp314 = _tmp309 + _tmp311 - _tmp313;
  const Scalar _tmp315 = _tmp308 * _tmp87;
  const Scalar _tmp316 = -_tmp107 * _tmp314 + _tmp130 * _tmp315 - _tmp131 * _tmp304;
  const Scalar _tmp317 = _tmp210 * _tmp308;
  const Scalar _tmp318 = _tmp102 * (-_tmp193 * _tmp305 + _tmp194 * _tmp308);
  const Scalar _tmp319 = _tmp100 * _tmp312;
  const Scalar _tmp320 = _tmp100 * _tmp310;
  const Scalar _tmp321 = _tmp318 - _tmp319 + _tmp320;
  const Scalar _tmp322 = _tmp104 * _tmp315 - _tmp107 * _tmp321 - _tmp200 * _tmp304;
  const Scalar _tmp323 = _tmp218 * _tmp305;
  const Scalar _tmp324 = _tmp121 * _tmp304;
  const Scalar _tmp325 =
      -_tmp117 * (-_tmp104 * _tmp317 + _tmp109 * _tmp322 + _tmp171 * _tmp200 + _tmp214 * _tmp321) -
      _tmp125 * (_tmp118 * _tmp171 - _tmp220 * _tmp305 + _tmp323 * _tmp58 - _tmp324 * _tmp58) -
      _tmp135 * (_tmp109 * _tmp316 - _tmp130 * _tmp317 + _tmp131 * _tmp171 + _tmp214 * _tmp314) -
      _tmp141 * (-_tmp137 * _tmp317 + _tmp138 * _tmp315 + _tmp171 * _tmp207 - _tmp211 * _tmp304);
  const Scalar _tmp326 = Scalar(1.4083112389913199) * _tmp325;
  const Scalar _tmp327 = _tmp190 * _tmp307;
  const Scalar _tmp328 = _tmp208 * (_tmp137 * _tmp327 - _tmp207 * _tmp306);
  const Scalar _tmp329 = _tmp147 * _tmp312;
  const Scalar _tmp330 = _tmp147 * _tmp310;
  const Scalar _tmp331 =
      -_tmp118 * _tmp306 - _tmp144 * _tmp329 + _tmp144 * _tmp330 + _tmp204 * _tmp305;
  const Scalar _tmp332 = _tmp130 * _tmp327 - _tmp131 * _tmp306 - _tmp157 * _tmp314 +
                         _tmp158 * _tmp309 + _tmp158 * _tmp311 - _tmp158 * _tmp313;
  const Scalar _tmp333 = _tmp104 * _tmp327 - _tmp157 * _tmp321 + _tmp158 * _tmp318 -
                         _tmp158 * _tmp319 + _tmp158 * _tmp320 - _tmp200 * _tmp306;
  const Scalar _tmp334 = _tmp167 * (_tmp168 * (_tmp151 * (-_tmp150 * _tmp331 - _tmp329 + _tmp330) +
                                               _tmp161 * (_tmp147 * _tmp309 + _tmp147 * _tmp311 -
                                                          _tmp147 * _tmp313 - _tmp150 * _tmp332) +
                                               _tmp163 * (_tmp147 * _tmp318 - _tmp147 * _tmp319 +
                                                          _tmp147 * _tmp320 - _tmp150 * _tmp333) -
                                               Scalar(1.0) * _tmp328) -
                                    _tmp224 * _tmp325);
  const Scalar _tmp335 = _tmp124 * _tmp323 - _tmp124 * _tmp324 - _tmp240 * _tmp304 +
                         _tmp242 * _tmp315 + _tmp244 * _tmp316 + _tmp245 * _tmp322;
  const Scalar _tmp336 =
      _tmp268 * (_tmp248 * (_tmp250 * _tmp332 + _tmp252 * _tmp333 + _tmp253 * _tmp331 + _tmp328) -
                 _tmp269 * _tmp335);
  const Scalar _tmp337 = -_tmp104 * _tmp116 * _tmp308 - _tmp130 * _tmp134 * _tmp308 -
                         _tmp241 * _tmp308 - _tmp273 * _tmp305 + _tmp274 * _tmp321 +
                         _tmp275 * _tmp314;
  const Scalar _tmp338 =
      _tmp296 *
      (_tmp278 * (_tmp124 * _tmp329 - _tmp124 * _tmp330 - _tmp280 * _tmp318 + _tmp280 * _tmp319 -
                  _tmp280 * _tmp320 - _tmp281 * _tmp309 - _tmp281 * _tmp311 + _tmp281 * _tmp313) -
       _tmp297 * _tmp337);

  // Output terms (1)
  Eigen::Matrix<Scalar, 4, 3> _res;

  _res(0, 0) = 0;
  _res(1, 0) =
      -Scalar(1.4083112389913199) * _tmp223 * _tmp238 -
      _tmp237 *
          (-_tmp223 * _tmp230 + _tmp225 * _tmp228 -
           _tmp235 * (-_tmp223 * _tmp236 + _tmp233 * (-_tmp223 * _tmp231 - _tmp225 * _tmp237)));
  _res(2, 0) =
      -_tmp246 * _tmp262 -
      _tmp257 *
          (-_tmp246 * _tmp272 -
           _tmp263 * (-_tmp246 * _tmp266 + _tmp259 * (-_tmp246 * _tmp267 - _tmp257 * _tmp270)) +
           _tmp270 * _tmp271);
  _res(3, 0) =
      -_tmp276 * _tmp291 -
      _tmp295 * (-_tmp276 * _tmp294 + _tmp298 * _tmp302 -
                 _tmp301 * (-_tmp276 * _tmp300 + _tmp279 * (-_tmp276 * _tmp284 - _tmp286 * _tmp299 -
                                                            _tmp295 * _tmp298)));
  _res(0, 1) = 0;
  _res(1, 1) =
      -_tmp237 *
          (_tmp228 * _tmp334 - _tmp230 * _tmp325 -
           _tmp235 * (_tmp233 * (-_tmp226 * _tmp326 - _tmp237 * _tmp334) - _tmp236 * _tmp325)) -
      _tmp238 * _tmp326;
  _res(2, 1) =
      -_tmp257 *
          (-_tmp263 * (_tmp259 * (-_tmp257 * _tmp336 - _tmp267 * _tmp335) - _tmp266 * _tmp335) +
           _tmp271 * _tmp336 - _tmp272 * _tmp335) -
      _tmp262 * _tmp335;
  _res(3, 1) =
      -_tmp291 * _tmp337 -
      _tmp295 * (-_tmp294 * _tmp337 -
                 _tmp301 * (_tmp279 * (-_tmp284 * _tmp337 - _tmp285 * _tmp299 - _tmp295 * _tmp338) -
                            _tmp300 * _tmp337) +
                 _tmp302 * _tmp338);
  _res(0, 2) = 0;
  _res(1, 2) = 0;
  _res(2, 2) = 0;
  _res(3, 2) = Scalar(-1.0);

  return _res;
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym

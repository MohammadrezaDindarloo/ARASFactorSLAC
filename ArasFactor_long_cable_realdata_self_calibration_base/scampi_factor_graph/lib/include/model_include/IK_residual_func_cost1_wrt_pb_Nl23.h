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
 * Symbolic function: IK_residual_func_cost1_wrt_pb_Nl23
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
Eigen::Matrix<Scalar, 4, 3> IkResidualFuncCost1WrtPbNl23(
    const Scalar fh1, const Scalar fv1, const sym::Rot3<Scalar>& DeltaRot,
    const Eigen::Matrix<Scalar, 3, 1>& position_vector, const Eigen::Matrix<Scalar, 4, 1>& encoder,
    const Eigen::Matrix<Scalar, 3, 1>& p_a, const Eigen::Matrix<Scalar, 3, 1>& p_b,
    const Eigen::Matrix<Scalar, 3, 1>& p_c, const Eigen::Matrix<Scalar, 3, 1>& p_d,
    const sym::Rot3<Scalar>& Rot_init, const Scalar epsilon) {
  // Total ops: 1166

  // Unused inputs
  (void)encoder;
  (void)epsilon;

  // Input arrays
  const Eigen::Matrix<Scalar, 4, 1>& _DeltaRot = DeltaRot.Data();
  const Eigen::Matrix<Scalar, 4, 1>& _Rot_init = Rot_init.Data();

  // Intermediate terms (368)
  const Scalar _tmp0 = _DeltaRot[0] * _Rot_init[2] + _DeltaRot[1] * _Rot_init[3] -
                       _DeltaRot[2] * _Rot_init[0] + _DeltaRot[3] * _Rot_init[1];
  const Scalar _tmp1 = _DeltaRot[0] * _Rot_init[3] - _DeltaRot[1] * _Rot_init[2] +
                       _DeltaRot[2] * _Rot_init[1] + _DeltaRot[3] * _Rot_init[0];
  const Scalar _tmp2 = 2 * _tmp0 * _tmp1;
  const Scalar _tmp3 = -_DeltaRot[0] * _Rot_init[1] + _DeltaRot[1] * _Rot_init[0] +
                       _DeltaRot[2] * _Rot_init[3] + _DeltaRot[3] * _Rot_init[2];
  const Scalar _tmp4 = -2 * _DeltaRot[0] * _Rot_init[0] - 2 * _DeltaRot[1] * _Rot_init[1] -
                       2 * _DeltaRot[2] * _Rot_init[2] + 2 * _DeltaRot[3] * _Rot_init[3];
  const Scalar _tmp5 = _tmp3 * _tmp4;
  const Scalar _tmp6 = Scalar(0.20999999999999999) * _tmp2 + Scalar(0.20999999999999999) * _tmp5;
  const Scalar _tmp7 = -_tmp6;
  const Scalar _tmp8 = 2 * _tmp3;
  const Scalar _tmp9 = _tmp0 * _tmp8;
  const Scalar _tmp10 = _tmp1 * _tmp4;
  const Scalar _tmp11 = -_tmp10 + _tmp9;
  const Scalar _tmp12 = -Scalar(0.010999999999999999) * _tmp11;
  const Scalar _tmp13 = -2 * std::pow(_tmp1, Scalar(2));
  const Scalar _tmp14 = 1 - 2 * std::pow(_tmp3, Scalar(2));
  const Scalar _tmp15 = Scalar(0.20999999999999999) * _tmp13 + Scalar(0.20999999999999999) * _tmp14;
  const Scalar _tmp16 = _tmp12 - _tmp15;
  const Scalar _tmp17 = _tmp16 + _tmp7;
  const Scalar _tmp18 = _tmp17 + position_vector(1, 0);
  const Scalar _tmp19 = _tmp18 - p_a(1, 0);
  const Scalar _tmp20 = -2 * std::pow(_tmp0, Scalar(2));
  const Scalar _tmp21 = Scalar(0.20999999999999999) * _tmp14 + Scalar(0.20999999999999999) * _tmp20;
  const Scalar _tmp22 = -_tmp21;
  const Scalar _tmp23 = _tmp1 * _tmp8;
  const Scalar _tmp24 = _tmp0 * _tmp4;
  const Scalar _tmp25 = _tmp23 + _tmp24;
  const Scalar _tmp26 = -Scalar(0.010999999999999999) * _tmp25;
  const Scalar _tmp27 = Scalar(0.20999999999999999) * _tmp2 - Scalar(0.20999999999999999) * _tmp5;
  const Scalar _tmp28 = _tmp26 - _tmp27;
  const Scalar _tmp29 = _tmp22 + _tmp28;
  const Scalar _tmp30 = _tmp29 + position_vector(0, 0);
  const Scalar _tmp31 = _tmp30 - p_a(0, 0);
  const Scalar _tmp32 = std::pow(Scalar(std::pow(_tmp19, Scalar(2)) + std::pow(_tmp31, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp33 = _tmp19 * _tmp32;
  const Scalar _tmp34 = _tmp26 + _tmp27;
  const Scalar _tmp35 = _tmp21 + _tmp34;
  const Scalar _tmp36 = _tmp35 + position_vector(0, 0);
  const Scalar _tmp37 = _tmp36 - p_c(0, 0);
  const Scalar _tmp38 = Scalar(1.0) / (_tmp37);
  const Scalar _tmp39 = _tmp12 + _tmp15;
  const Scalar _tmp40 = _tmp39 + _tmp6;
  const Scalar _tmp41 = _tmp40 + position_vector(1, 0);
  const Scalar _tmp42 = _tmp41 - p_c(1, 0);
  const Scalar _tmp43 = _tmp38 * _tmp42;
  const Scalar _tmp44 = _tmp31 * _tmp32;
  const Scalar _tmp45 = -_tmp33 + _tmp43 * _tmp44;
  const Scalar _tmp46 = _tmp21 + _tmp28;
  const Scalar _tmp47 = _tmp46 + position_vector(0, 0);
  const Scalar _tmp48 = _tmp47 - p_b(0, 0);
  const Scalar _tmp49 = std::pow(_tmp48, Scalar(2));
  const Scalar _tmp50 = _tmp16 + _tmp6;
  const Scalar _tmp51 = _tmp50 + position_vector(1, 0);
  const Scalar _tmp52 = _tmp51 - p_b(1, 0);
  const Scalar _tmp53 = std::pow(_tmp52, Scalar(2));
  const Scalar _tmp54 = _tmp49 + _tmp53;
  const Scalar _tmp55 = std::pow(_tmp54, Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp56 = _tmp43 * _tmp55;
  const Scalar _tmp57 = _tmp48 * _tmp56 - _tmp52 * _tmp55;
  const Scalar _tmp58 = std::pow(_tmp57, Scalar(-2));
  const Scalar _tmp59 = _tmp45 * _tmp58;
  const Scalar _tmp60 = Scalar(0.20999999999999999) * _tmp10 + Scalar(0.20999999999999999) * _tmp9;
  const Scalar _tmp61 = -_tmp60;
  const Scalar _tmp62 = -Scalar(0.010999999999999999) * _tmp13 -
                        Scalar(0.010999999999999999) * _tmp20 + Scalar(-0.010999999999999999);
  const Scalar _tmp63 = Scalar(0.20999999999999999) * _tmp23 - Scalar(0.20999999999999999) * _tmp24;
  const Scalar _tmp64 = _tmp62 + _tmp63;
  const Scalar _tmp65 = _tmp61 + _tmp64;
  const Scalar _tmp66 = _tmp55 * _tmp65;
  const Scalar _tmp67 = _tmp60 + _tmp64;
  const Scalar _tmp68 = _tmp55 * _tmp67;
  const Scalar _tmp69 = _tmp48 * _tmp68;
  const Scalar _tmp70 = -_tmp48 * _tmp66 + _tmp69;
  const Scalar _tmp71 = std::pow(_tmp54, Scalar(Scalar(-3) / Scalar(2)));
  const Scalar _tmp72 = _tmp48 * _tmp52 * _tmp71;
  const Scalar _tmp73 = _tmp49 * _tmp71;
  const Scalar _tmp74 = _tmp43 * _tmp73 - _tmp56 - _tmp72;
  const Scalar _tmp75 = _tmp70 * _tmp74;
  const Scalar _tmp76 = _tmp67 * _tmp73;
  const Scalar _tmp77 = -_tmp65 * _tmp73 + _tmp66 - _tmp68 + _tmp76;
  const Scalar _tmp78 = Scalar(1.0) / (_tmp57);
  const Scalar _tmp79 = _tmp45 * _tmp78;
  const Scalar _tmp80 = _tmp65 * _tmp72;
  const Scalar _tmp81 = _tmp43 * _tmp68 - _tmp43 * _tmp76 + _tmp80;
  const Scalar _tmp82 = -_tmp43 * _tmp69 + _tmp52 * _tmp66;
  const Scalar _tmp83 = _tmp74 * _tmp82;
  const Scalar _tmp84 = _tmp59 * _tmp83 - _tmp79 * _tmp81;
  const Scalar _tmp85 = Scalar(1.0) * _tmp40;
  const Scalar _tmp86 = -_tmp85;
  const Scalar _tmp87 = Scalar(1.0) / (_tmp50 + _tmp86);
  const Scalar _tmp88 = Scalar(1.0) * _tmp35;
  const Scalar _tmp89 = -_tmp46 + _tmp88;
  const Scalar _tmp90 = _tmp87 * _tmp89;
  const Scalar _tmp91 = _tmp59 * _tmp75 - _tmp77 * _tmp79 - _tmp84 * _tmp90;
  const Scalar _tmp92 = _tmp62 - _tmp63;
  const Scalar _tmp93 = _tmp61 + _tmp92;
  const Scalar _tmp94 = _tmp44 * _tmp67;
  const Scalar _tmp95 = _tmp33 * _tmp93 - _tmp43 * _tmp94 - _tmp79 * _tmp82;
  const Scalar _tmp96 = -_tmp44 * _tmp93 - _tmp70 * _tmp79 - _tmp90 * _tmp95 + _tmp94;
  const Scalar _tmp97 = std::pow(_tmp96, Scalar(-2));
  const Scalar _tmp98 = _tmp91 * _tmp97;
  const Scalar _tmp99 = _tmp45 * _tmp98;
  const Scalar _tmp100 = _tmp85 * _tmp90 + _tmp88;
  const Scalar _tmp101 = 0;
  const Scalar _tmp102 = _tmp55 * _tmp78;
  const Scalar _tmp103 = _tmp102 * _tmp48;
  const Scalar _tmp104 = _tmp101 * _tmp103;
  const Scalar _tmp105 = _tmp73 * _tmp78;
  const Scalar _tmp106 = Scalar(1.0) / (_tmp96);
  const Scalar _tmp107 = _tmp106 * _tmp45;
  const Scalar _tmp108 = _tmp101 * _tmp107;
  const Scalar _tmp109 = _tmp101 * _tmp98;
  const Scalar _tmp110 = _tmp58 * _tmp74;
  const Scalar _tmp111 = _tmp48 * _tmp55;
  const Scalar _tmp112 = _tmp108 * _tmp111;
  const Scalar _tmp113 = Scalar(333.54000000000002) - fv1;
  const Scalar _tmp114 =
      std::sqrt(Scalar(std::pow(_tmp37, Scalar(2)) + std::pow(_tmp42, Scalar(2))));
  const Scalar _tmp115 = _tmp114 * _tmp38;
  const Scalar _tmp116 = _tmp113 * _tmp115;
  const Scalar _tmp117 = Scalar(1.0) / (_tmp114);
  const Scalar _tmp118 = _tmp115 * (_tmp117 * _tmp35 * _tmp42 - _tmp117 * _tmp37 * _tmp40);
  const Scalar _tmp119 = _tmp118 * _tmp55;
  const Scalar _tmp120 = _tmp50 * _tmp55;
  const Scalar _tmp121 = _tmp118 * _tmp73 - _tmp119 - _tmp120 - _tmp46 * _tmp72 + _tmp50 * _tmp73;
  const Scalar _tmp122 = _tmp46 * _tmp55;
  const Scalar _tmp123 = _tmp119 * _tmp48 + _tmp120 * _tmp48 - _tmp122 * _tmp52;
  const Scalar _tmp124 = _tmp123 * _tmp74;
  const Scalar _tmp125 = -_tmp121 * _tmp79 + _tmp124 * _tmp59;
  const Scalar _tmp126 = _tmp118 * _tmp44 - _tmp123 * _tmp79 + _tmp17 * _tmp44 - _tmp29 * _tmp33;
  const Scalar _tmp127 = std::pow(_tmp126, Scalar(-2));
  const Scalar _tmp128 = _tmp125 * _tmp127;
  const Scalar _tmp129 = Scalar(1.0) * _tmp103 * _tmp45;
  const Scalar _tmp130 = Scalar(1.0) / (_tmp126);
  const Scalar _tmp131 = Scalar(1.0) * _tmp130;
  const Scalar _tmp132 = _tmp131 * _tmp45;
  const Scalar _tmp133 = _tmp131 * _tmp79;
  const Scalar _tmp134 = _tmp59 * _tmp74;
  const Scalar _tmp135 = _tmp111 * _tmp131;
  const Scalar _tmp136 = Scalar(1.0) * _tmp44;
  const Scalar _tmp137 = _tmp22 + _tmp34;
  const Scalar _tmp138 = _tmp137 - p_d(0, 0) + position_vector(0, 0);
  const Scalar _tmp139 = _tmp39 + _tmp7;
  const Scalar _tmp140 = _tmp139 - p_d(1, 0) + position_vector(1, 0);
  const Scalar _tmp141 =
      std::pow(Scalar(std::pow(_tmp138, Scalar(2)) + std::pow(_tmp140, Scalar(2))),
               Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp142 = _tmp140 * _tmp141;
  const Scalar _tmp143 = _tmp138 * _tmp141;
  const Scalar _tmp144 = fh1 * (_tmp137 * _tmp142 - _tmp139 * _tmp143);
  const Scalar _tmp145 = _tmp115 * _tmp144;
  const Scalar _tmp146 = _tmp43 * _tmp67;
  const Scalar _tmp147 = _tmp43 * _tmp78;
  const Scalar _tmp148 = _tmp146 + _tmp147 * _tmp82;
  const Scalar _tmp149 = _tmp147 * _tmp70 - _tmp148 * _tmp90 - _tmp67;
  const Scalar _tmp150 = _tmp106 * _tmp126;
  const Scalar _tmp151 = -_tmp118 + _tmp123 * _tmp147 - _tmp149 * _tmp150;
  const Scalar _tmp152 = _tmp130 * _tmp96;
  const Scalar _tmp153 = _tmp151 * _tmp152;
  const Scalar _tmp154 = _tmp149 + _tmp153;
  const Scalar _tmp155 = -_tmp107 * _tmp154 - _tmp43;
  const Scalar _tmp156 = _tmp130 * _tmp151;
  const Scalar _tmp157 = _tmp156 * _tmp91;
  const Scalar _tmp158 = _tmp106 * _tmp149;
  const Scalar _tmp159 = _tmp126 * _tmp98;
  const Scalar _tmp160 = _tmp43 * _tmp58;
  const Scalar _tmp161 = _tmp147 * _tmp81 - _tmp160 * _tmp83;
  const Scalar _tmp162 = _tmp147 * _tmp77 - _tmp160 * _tmp75 - _tmp161 * _tmp90;
  const Scalar _tmp163 = _tmp152 * (_tmp121 * _tmp147 - _tmp124 * _tmp160 - _tmp125 * _tmp158 +
                                    _tmp149 * _tmp159 - _tmp150 * _tmp162);
  const Scalar _tmp164 = _tmp151 * _tmp96;
  const Scalar _tmp165 = _tmp128 * _tmp164;
  const Scalar _tmp166 = _tmp157 + _tmp162 + _tmp163 - _tmp165;
  const Scalar _tmp167 = _tmp106 * _tmp44;
  const Scalar _tmp168 = _tmp102 * _tmp155;
  const Scalar _tmp169 = _tmp111 * _tmp155;
  const Scalar _tmp170 = -_tmp107 * _tmp166 + _tmp154 * _tmp99;
  const Scalar _tmp171 = _tmp154 * _tmp44;
  const Scalar _tmp172 = _tmp143 * fh1;
  const Scalar _tmp173 = _tmp115 * _tmp172;
  const Scalar _tmp174 = Scalar(1.0) * _tmp78;
  const Scalar _tmp175 = Scalar(1.0) * _tmp87;
  const Scalar _tmp176 = _tmp175 * _tmp89;
  const Scalar _tmp177 = _tmp176 * _tmp78;
  const Scalar _tmp178 = -_tmp174 * _tmp70 + _tmp177 * _tmp82;
  const Scalar _tmp179 = _tmp106 * _tmp178;
  const Scalar _tmp180 = -_tmp123 * _tmp174 - _tmp126 * _tmp179;
  const Scalar _tmp181 = _tmp152 * _tmp180;
  const Scalar _tmp182 = _tmp178 + _tmp181;
  const Scalar _tmp183 = -_tmp107 * _tmp182 + Scalar(1.0);
  const Scalar _tmp184 = _tmp180 * _tmp96;
  const Scalar _tmp185 = _tmp128 * _tmp184;
  const Scalar _tmp186 = _tmp176 * _tmp58;
  const Scalar _tmp187 = Scalar(1.0) * _tmp58;
  const Scalar _tmp188 = -_tmp174 * _tmp77 + _tmp177 * _tmp81 - _tmp186 * _tmp83 + _tmp187 * _tmp75;
  const Scalar _tmp189 = _tmp152 * (-_tmp121 * _tmp174 + _tmp124 * _tmp187 - _tmp125 * _tmp179 -
                                    _tmp150 * _tmp188 + _tmp159 * _tmp178);
  const Scalar _tmp190 = _tmp130 * _tmp180;
  const Scalar _tmp191 = _tmp190 * _tmp91;
  const Scalar _tmp192 = -_tmp185 + _tmp188 + _tmp189 + _tmp191;
  const Scalar _tmp193 = -_tmp107 * _tmp192 + _tmp182 * _tmp99;
  const Scalar _tmp194 = _tmp110 * _tmp183;
  const Scalar _tmp195 = _tmp182 * _tmp44;
  const Scalar _tmp196 = _tmp102 * _tmp183;
  const Scalar _tmp197 = _tmp142 * fh1;
  const Scalar _tmp198 = _tmp115 * _tmp197;
  const Scalar _tmp199 = -_tmp116 * (_tmp102 * _tmp108 + _tmp104 * _tmp99 - _tmp105 * _tmp108 -
                                     _tmp109 * _tmp44 + _tmp110 * _tmp112) -
                         _tmp145 * (_tmp102 * _tmp132 + _tmp128 * _tmp129 - _tmp128 * _tmp136 -
                                    _tmp133 * _tmp73 + _tmp134 * _tmp135) -
                         _tmp173 * (_tmp103 * _tmp170 + _tmp105 * _tmp155 - _tmp110 * _tmp169 +
                                    _tmp166 * _tmp167 - _tmp168 - _tmp171 * _tmp98) -
                         _tmp198 * (_tmp103 * _tmp193 + _tmp105 * _tmp183 - _tmp111 * _tmp194 +
                                    _tmp167 * _tmp192 - _tmp195 * _tmp98 - _tmp196);
  const Scalar _tmp200 = _tmp101 * _tmp106;
  const Scalar _tmp201 = -_tmp116 * (-_tmp103 * _tmp108 + _tmp200 * _tmp44) -
                         _tmp145 * (-_tmp103 * _tmp132 + _tmp131 * _tmp44) -
                         _tmp173 * (_tmp154 * _tmp167 + _tmp168 * _tmp48 + Scalar(1.0)) -
                         _tmp198 * (_tmp167 * _tmp182 + _tmp196 * _tmp48);
  const Scalar _tmp202 = std::pow(_tmp201, Scalar(-2));
  const Scalar _tmp203 = Scalar(0.71007031138673404) * _tmp202;
  const Scalar _tmp204 = _tmp199 * _tmp203;
  const Scalar _tmp205 = _tmp17 + _tmp86;
  const Scalar _tmp206 = _tmp205 * _tmp90;
  const Scalar _tmp207 = Scalar(1.0) / (-_tmp206 - _tmp29 + _tmp88);
  const Scalar _tmp208 = _tmp100 * _tmp207;
  const Scalar _tmp209 = _tmp106 * _tmp95;
  const Scalar _tmp210 = -_tmp101 * _tmp209 - _tmp205 * _tmp208 + _tmp86;
  const Scalar _tmp211 = Scalar(1.0) * _tmp207;
  const Scalar _tmp212 = _tmp205 * _tmp87;
  const Scalar _tmp213 = fh1 * (_tmp60 + _tmp92);
  const Scalar _tmp214 = _tmp137 * fv1 + _tmp143 * _tmp213 + Scalar(40.024799999999999) * _tmp25;
  const Scalar _tmp215 = _tmp205 * _tmp207;
  const Scalar _tmp216 = _tmp148 + _tmp153 * _tmp215 - _tmp154 * _tmp209;
  const Scalar _tmp217 = Scalar(1.0) * _tmp172;
  const Scalar _tmp218 = -_tmp174 * _tmp82 + _tmp181 * _tmp215 - _tmp182 * _tmp209;
  const Scalar _tmp219 = Scalar(1.0) * _tmp197;
  const Scalar _tmp220 = -Scalar(40.024799999999999) * _tmp11 - _tmp139 * fv1 - _tmp142 * _tmp213;
  const Scalar _tmp221 = _tmp206 * _tmp211 + Scalar(1.0);
  const Scalar _tmp222 = _tmp211 * _tmp90;
  const Scalar _tmp223 = _tmp152 * _tmp211;
  const Scalar _tmp224 = -_tmp131 * _tmp95 + _tmp205 * _tmp223;
  const Scalar _tmp225 = Scalar(1.0) * _tmp144;
  const Scalar _tmp226 =
      Scalar(1.0) * _tmp113 * (-_tmp100 * _tmp211 - _tmp175 * _tmp210 + Scalar(1.0)) +
      Scalar(1.0) * _tmp214 * (_tmp211 * _tmp212 - _tmp211) +
      _tmp217 * (_tmp153 * _tmp211 - _tmp175 * _tmp216) +
      _tmp219 * (-_tmp175 * _tmp218 + _tmp181 * _tmp211) +
      Scalar(1.0) * _tmp220 * (-_tmp175 * _tmp221 + _tmp222) +
      _tmp225 * (-_tmp175 * _tmp224 + _tmp223);
  const Scalar _tmp227 = Scalar(1.0) / (_tmp201);
  const Scalar _tmp228 = std::asinh(_tmp226 * _tmp227);
  const Scalar _tmp229 = Scalar(1.4083112389913199) * _tmp228;
  const Scalar _tmp230 = Scalar(1.4083112389913199) * _tmp201;
  const Scalar _tmp231 =
      std::pow(Scalar(_tmp202 * std::pow(_tmp226, Scalar(2)) + 1), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp232 = _tmp211 * _tmp96;
  const Scalar _tmp233 = _tmp128 * _tmp232;
  const Scalar _tmp234 = _tmp130 * _tmp211;
  const Scalar _tmp235 = _tmp234 * _tmp91;
  const Scalar _tmp236 = Scalar(1.0) * _tmp95;
  const Scalar _tmp237 =
      _tmp128 * _tmp236 - _tmp131 * _tmp84 - _tmp205 * _tmp233 + _tmp205 * _tmp235;
  const Scalar _tmp238 = _tmp154 * _tmp95;
  const Scalar _tmp239 = _tmp106 * _tmp84;
  const Scalar _tmp240 = -_tmp154 * _tmp239 + _tmp157 * _tmp215 + _tmp161 + _tmp163 * _tmp215 -
                         _tmp165 * _tmp215 - _tmp166 * _tmp209 + _tmp238 * _tmp98;
  const Scalar _tmp241 = _tmp182 * _tmp95;
  const Scalar _tmp242 = -_tmp174 * _tmp81 - _tmp182 * _tmp239 - _tmp185 * _tmp215 +
                         _tmp187 * _tmp83 + _tmp189 * _tmp215 + _tmp191 * _tmp215 -
                         _tmp192 * _tmp209 + _tmp241 * _tmp98;
  const Scalar _tmp243 = -_tmp101 * _tmp239 + _tmp109 * _tmp95;
  const Scalar _tmp244 = _tmp113 * _tmp175;
  const Scalar _tmp245 = _tmp202 * _tmp226;
  const Scalar _tmp246 = _tmp231 * (-_tmp199 * _tmp245 +
                                    _tmp227 * (_tmp217 * (_tmp157 * _tmp211 + _tmp163 * _tmp211 -
                                                          _tmp165 * _tmp211 - _tmp175 * _tmp240) +
                                               _tmp219 * (-_tmp175 * _tmp242 - _tmp185 * _tmp211 +
                                                          _tmp189 * _tmp211 + _tmp191 * _tmp211) +
                                               _tmp225 * (-_tmp175 * _tmp237 - _tmp233 + _tmp235) -
                                               _tmp243 * _tmp244));
  const Scalar _tmp247 = Scalar(0.71007031138673404) * _tmp227;
  const Scalar _tmp248 =
      -_tmp201 * _tmp229 - std::sqrt(Scalar(std::pow(Scalar(-_tmp36 + p_c(0, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp41 + p_c(1, 0)), Scalar(2))));
  const Scalar _tmp249 = _tmp247 * _tmp248;
  const Scalar _tmp250 = std::sinh(_tmp249);
  const Scalar _tmp251 = Scalar(1.0) * _tmp228;
  const Scalar _tmp252 = Scalar(1.0) * std::sinh(_tmp251);
  const Scalar _tmp253 = Scalar(1.4083112389913199) * _tmp247 * p_c(2, 0) -
                         Scalar(1.4083112389913199) * std::cosh(_tmp249) +
                         Scalar(1.4083112389913199) * std::cosh(_tmp251);
  const Scalar _tmp254 = _tmp108 * _tmp113;
  const Scalar _tmp255 = _tmp172 * _tmp78;
  const Scalar _tmp256 = _tmp131 * _tmp144;
  const Scalar _tmp257 = _tmp128 * _tmp225;
  const Scalar _tmp258 = _tmp155 * _tmp172;
  const Scalar _tmp259 = _tmp197 * _tmp78;
  const Scalar _tmp260 = _tmp109 * _tmp113;
  const Scalar _tmp261 = _tmp110 * _tmp254 - _tmp110 * _tmp258 + _tmp134 * _tmp256 +
                         _tmp170 * _tmp255 + _tmp193 * _tmp259 - _tmp194 * _tmp197 +
                         _tmp257 * _tmp79 + _tmp260 * _tmp79;
  const Scalar _tmp262 = _tmp113 * _tmp87;
  const Scalar _tmp263 = _tmp197 * _tmp87;
  const Scalar _tmp264 = _tmp144 * _tmp87;
  const Scalar _tmp265 = _tmp172 * _tmp87;
  const Scalar _tmp266 = _tmp211 * _tmp214;
  const Scalar _tmp267 = _tmp210 * _tmp262 - _tmp212 * _tmp266 + _tmp216 * _tmp265 +
                         _tmp218 * _tmp263 + _tmp220 * _tmp221 * _tmp87 + _tmp224 * _tmp264;
  const Scalar _tmp268 =
      _tmp155 * _tmp255 + _tmp183 * _tmp259 - _tmp254 * _tmp78 - _tmp256 * _tmp79;
  const Scalar _tmp269 = Scalar(1.0) / (_tmp268);
  const Scalar _tmp270 = std::asinh(_tmp267 * _tmp269);
  const Scalar _tmp271 = Scalar(1.0) * _tmp270;
  const Scalar _tmp272 = Scalar(0.71007031138673404) * _tmp269;
  const Scalar _tmp273 = -_tmp51 + p_b(1, 0);
  const Scalar _tmp274 = -_tmp47 + p_b(0, 0);
  const Scalar _tmp275 =
      std::sqrt(Scalar(std::pow(_tmp273, Scalar(2)) + std::pow(_tmp274, Scalar(2))));
  const Scalar _tmp276 = Scalar(1.4083112389913199) * _tmp270;
  const Scalar _tmp277 = -_tmp268 * _tmp276 - _tmp275;
  const Scalar _tmp278 = _tmp272 * _tmp277;
  const Scalar _tmp279 = Scalar(1.4083112389913199) * _tmp272 * p_b(2, 0) +
                         Scalar(1.4083112389913199) * std::cosh(_tmp271) -
                         Scalar(1.4083112389913199) * std::cosh(_tmp278);
  const Scalar _tmp280 = std::pow(_tmp268, Scalar(-2));
  const Scalar _tmp281 = _tmp267 * _tmp280;
  const Scalar _tmp282 = -_tmp261 * _tmp281 + _tmp269 * (_tmp237 * _tmp264 + _tmp240 * _tmp265 +
                                                         _tmp242 * _tmp263 + _tmp243 * _tmp262);
  const Scalar _tmp283 =
      std::pow(Scalar(std::pow(_tmp267, Scalar(2)) * _tmp280 + 1), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp284 = Scalar(1.0) * _tmp283 * std::sinh(_tmp271);
  const Scalar _tmp285 = Scalar(0.71007031138673404) * _tmp280;
  const Scalar _tmp286 = _tmp261 * _tmp285;
  const Scalar _tmp287 = Scalar(1.4083112389913199) * _tmp268;
  const Scalar _tmp288 = _tmp283 * _tmp287;
  const Scalar _tmp289 = Scalar(1.0) / (_tmp275);
  const Scalar _tmp290 = std::sinh(_tmp278);
  const Scalar _tmp291 = _tmp197 * _tmp207;
  const Scalar _tmp292 = _tmp172 * _tmp207;
  const Scalar _tmp293 = _tmp106 * _tmp197;
  const Scalar _tmp294 = _tmp106 * _tmp172;
  const Scalar _tmp295 = _tmp113 * _tmp200 + _tmp154 * _tmp294 + _tmp182 * _tmp293 + _tmp256;
  const Scalar _tmp296 = Scalar(1.0) / (_tmp295);
  const Scalar _tmp297 = _tmp182 * _tmp197;
  const Scalar _tmp298 = _tmp154 * _tmp172;
  const Scalar _tmp299 = _tmp166 * _tmp294 + _tmp192 * _tmp293 - _tmp257 - _tmp260 -
                         _tmp297 * _tmp98 - _tmp298 * _tmp98;
  const Scalar _tmp300 = std::pow(_tmp295, Scalar(-2));
  const Scalar _tmp301 = _tmp113 * _tmp208 - _tmp144 * _tmp223 - _tmp153 * _tmp292 -
                         _tmp181 * _tmp291 - _tmp220 * _tmp222 + _tmp266;
  const Scalar _tmp302 = _tmp300 * _tmp301;
  const Scalar _tmp303 =
      std::pow(Scalar(_tmp300 * std::pow(_tmp301, Scalar(2)) + 1), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp304 =
      _tmp303 *
      (_tmp296 * (_tmp144 * _tmp233 - _tmp144 * _tmp235 - _tmp157 * _tmp292 - _tmp163 * _tmp292 +
                  _tmp165 * _tmp292 + _tmp185 * _tmp291 - _tmp189 * _tmp291 - _tmp191 * _tmp291) -
       _tmp299 * _tmp302);
  const Scalar _tmp305 = std::asinh(_tmp296 * _tmp301);
  const Scalar _tmp306 = Scalar(1.0) * _tmp305;
  const Scalar _tmp307 = Scalar(1.0) * std::sinh(_tmp306);
  const Scalar _tmp308 = Scalar(1.4083112389913199) * _tmp305;
  const Scalar _tmp309 =
      -_tmp295 * _tmp308 - std::sqrt(Scalar(std::pow(Scalar(-_tmp18 + p_a(1, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp30 + p_a(0, 0)), Scalar(2))));
  const Scalar _tmp310 = Scalar(0.71007031138673404) * _tmp296;
  const Scalar _tmp311 = _tmp309 * _tmp310;
  const Scalar _tmp312 = std::sinh(_tmp311);
  const Scalar _tmp313 = Scalar(1.4083112389913199) * _tmp295;
  const Scalar _tmp314 = Scalar(0.71007031138673404) * _tmp300;
  const Scalar _tmp315 = _tmp309 * _tmp314;
  const Scalar _tmp316 = _tmp314 * p_a(2, 0);
  const Scalar _tmp317 = Scalar(1.4083112389913199) * _tmp310 * p_a(2, 0) +
                         Scalar(1.4083112389913199) * std::cosh(_tmp306) -
                         Scalar(1.4083112389913199) * std::cosh(_tmp311);
  const Scalar _tmp318 = _tmp53 * _tmp71;
  const Scalar _tmp319 = -_tmp318 + _tmp43 * _tmp72 + _tmp55;
  const Scalar _tmp320 = _tmp123 * _tmp319;
  const Scalar _tmp321 = _tmp78 * (_tmp118 * _tmp72 + _tmp122 - _tmp318 * _tmp46 + _tmp50 * _tmp72);
  const Scalar _tmp322 = _tmp320 * _tmp59 - _tmp321 * _tmp45;
  const Scalar _tmp323 = _tmp127 * _tmp322;
  const Scalar _tmp324 = _tmp319 * _tmp59;
  const Scalar _tmp325 = _tmp319 * _tmp70;
  const Scalar _tmp326 = _tmp319 * _tmp82;
  const Scalar _tmp327 = -_tmp146 * _tmp72 + _tmp318 * _tmp65 - _tmp66;
  const Scalar _tmp328 = _tmp326 * _tmp59 - _tmp327 * _tmp79;
  const Scalar _tmp329 = _tmp67 * _tmp72 - _tmp80;
  const Scalar _tmp330 = _tmp325 * _tmp59 - _tmp328 * _tmp90 - _tmp329 * _tmp79;
  const Scalar _tmp331 = _tmp330 * _tmp97;
  const Scalar _tmp332 = _tmp72 * _tmp78;
  const Scalar _tmp333 = _tmp319 * _tmp58;
  const Scalar _tmp334 = _tmp164 * _tmp323;
  const Scalar _tmp335 = _tmp156 * _tmp330;
  const Scalar _tmp336 = _tmp126 * _tmp331;
  const Scalar _tmp337 = _tmp147 * _tmp327 - _tmp160 * _tmp326;
  const Scalar _tmp338 = _tmp147 * _tmp329 - _tmp160 * _tmp325 - _tmp337 * _tmp90;
  const Scalar _tmp339 = _tmp152 * (_tmp149 * _tmp336 - _tmp150 * _tmp338 - _tmp158 * _tmp322 -
                                    _tmp160 * _tmp320 + _tmp321 * _tmp43);
  const Scalar _tmp340 = -_tmp334 + _tmp335 + _tmp338 + _tmp339;
  const Scalar _tmp341 = _tmp331 * _tmp45;
  const Scalar _tmp342 = -_tmp107 * _tmp340 + _tmp154 * _tmp341;
  const Scalar _tmp343 = _tmp101 * _tmp331;
  const Scalar _tmp344 =
      -_tmp174 * _tmp329 + _tmp177 * _tmp327 - _tmp186 * _tmp326 + _tmp187 * _tmp325;
  const Scalar _tmp345 = _tmp152 * (-_tmp150 * _tmp344 + _tmp178 * _tmp336 - _tmp179 * _tmp322 +
                                    _tmp187 * _tmp320 - Scalar(1.0) * _tmp321);
  const Scalar _tmp346 = _tmp184 * _tmp323;
  const Scalar _tmp347 = _tmp190 * _tmp330;
  const Scalar _tmp348 = _tmp344 + _tmp345 - _tmp346 + _tmp347;
  const Scalar _tmp349 = _tmp183 * _tmp333;
  const Scalar _tmp350 = -_tmp107 * _tmp348 + _tmp182 * _tmp341;
  const Scalar _tmp351 =
      -_tmp116 * (_tmp104 * _tmp341 - _tmp108 * _tmp332 + _tmp112 * _tmp333 - _tmp343 * _tmp44) -
      _tmp145 * (_tmp129 * _tmp323 - _tmp133 * _tmp72 + _tmp135 * _tmp324 - _tmp136 * _tmp323) -
      _tmp173 * (_tmp103 * _tmp342 + _tmp155 * _tmp332 + _tmp167 * _tmp340 - _tmp169 * _tmp333 -
                 _tmp171 * _tmp331) -
      _tmp198 * (_tmp103 * _tmp350 - _tmp111 * _tmp349 + _tmp167 * _tmp348 + _tmp183 * _tmp332 -
                 _tmp195 * _tmp331);
  const Scalar _tmp352 = _tmp232 * _tmp323;
  const Scalar _tmp353 = _tmp234 * _tmp330;
  const Scalar _tmp354 =
      -_tmp131 * _tmp328 - _tmp205 * _tmp352 + _tmp205 * _tmp353 + _tmp236 * _tmp323;
  const Scalar _tmp355 = _tmp106 * _tmp328;
  const Scalar _tmp356 = -_tmp154 * _tmp355 - _tmp209 * _tmp340 - _tmp215 * _tmp334 +
                         _tmp215 * _tmp335 + _tmp215 * _tmp339 + _tmp238 * _tmp331 + _tmp337;
  const Scalar _tmp357 = -_tmp174 * _tmp327 - _tmp182 * _tmp355 + _tmp187 * _tmp326 -
                         _tmp209 * _tmp348 + _tmp215 * _tmp345 - _tmp215 * _tmp346 +
                         _tmp215 * _tmp347 + _tmp241 * _tmp331;
  const Scalar _tmp358 = -_tmp101 * _tmp355 + _tmp343 * _tmp95;
  const Scalar _tmp359 = _tmp231 * (_tmp227 * (_tmp217 * (-_tmp175 * _tmp356 - _tmp211 * _tmp334 +
                                                          _tmp211 * _tmp335 + _tmp211 * _tmp339) +
                                               _tmp219 * (-_tmp175 * _tmp357 + _tmp211 * _tmp345 -
                                                          _tmp211 * _tmp346 + _tmp211 * _tmp347) +
                                               _tmp225 * (-_tmp175 * _tmp354 - _tmp352 + _tmp353) -
                                               _tmp244 * _tmp358) -
                                    _tmp245 * _tmp351);
  const Scalar _tmp360 = _tmp203 * _tmp351;
  const Scalar _tmp361 = _tmp225 * _tmp323;
  const Scalar _tmp362 = _tmp113 * _tmp343;
  const Scalar _tmp363 = -_tmp197 * _tmp349 + _tmp254 * _tmp333 + _tmp255 * _tmp342 +
                         _tmp256 * _tmp324 - _tmp258 * _tmp333 + _tmp259 * _tmp350 +
                         _tmp361 * _tmp79 + _tmp362 * _tmp79;
  const Scalar _tmp364 = _tmp285 * _tmp363;
  const Scalar _tmp365 =
      _tmp269 * (_tmp262 * _tmp358 + _tmp263 * _tmp357 + _tmp264 * _tmp354 + _tmp265 * _tmp356) -
      _tmp281 * _tmp363;
  const Scalar _tmp366 = _tmp293 * _tmp348 + _tmp294 * _tmp340 - _tmp297 * _tmp331 -
                         _tmp298 * _tmp331 - _tmp361 - _tmp362;
  const Scalar _tmp367 =
      _tmp303 *
      (_tmp296 * (_tmp144 * _tmp352 - _tmp144 * _tmp353 - _tmp291 * _tmp345 + _tmp291 * _tmp346 -
                  _tmp291 * _tmp347 + _tmp292 * _tmp334 - _tmp292 * _tmp335 - _tmp292 * _tmp339) -
       _tmp302 * _tmp366);

  // Output terms (1)
  Eigen::Matrix<Scalar, 4, 3> _res;

  _res(0, 0) = 0;
  _res(1, 0) =
      -_tmp199 * _tmp253 -
      _tmp230 *
          (-_tmp204 * p_c(2, 0) + _tmp246 * _tmp252 -
           _tmp250 * (-_tmp204 * _tmp248 + _tmp247 * (-_tmp199 * _tmp229 - _tmp230 * _tmp246)));
  _res(2, 0) =
      -_tmp261 * _tmp279 -
      _tmp287 * (_tmp282 * _tmp284 - _tmp286 * p_b(2, 0) -
                 _tmp290 * (_tmp272 * (-_tmp261 * _tmp276 - _tmp274 * _tmp289 - _tmp282 * _tmp288) -
                            _tmp277 * _tmp286));
  _res(3, 0) =
      -_tmp299 * _tmp317 -
      _tmp313 *
          (-_tmp299 * _tmp316 + _tmp304 * _tmp307 -
           _tmp312 * (-_tmp299 * _tmp315 + _tmp310 * (-_tmp299 * _tmp308 - _tmp304 * _tmp313)));
  _res(0, 1) = 0;
  _res(1, 1) =
      -_tmp230 *
          (-_tmp250 * (_tmp247 * (-_tmp229 * _tmp351 - _tmp230 * _tmp359) - _tmp248 * _tmp360) +
           _tmp252 * _tmp359 - _tmp360 * p_c(2, 0)) -
      _tmp253 * _tmp351;
  _res(2, 1) =
      -_tmp279 * _tmp363 -
      _tmp287 * (_tmp284 * _tmp365 -
                 _tmp290 * (_tmp272 * (-_tmp273 * _tmp289 - _tmp276 * _tmp363 - _tmp288 * _tmp365) -
                            _tmp277 * _tmp364) -
                 _tmp364 * p_b(2, 0));
  _res(3, 1) =
      -_tmp313 *
          (_tmp307 * _tmp367 -
           _tmp312 * (_tmp310 * (-_tmp308 * _tmp366 - _tmp313 * _tmp367) - _tmp315 * _tmp366) -
           _tmp316 * _tmp366) -
      _tmp317 * _tmp366;
  _res(0, 2) = 0;
  _res(1, 2) = 0;
  _res(2, 2) = Scalar(-1.0);
  _res(3, 2) = 0;

  return _res;
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym

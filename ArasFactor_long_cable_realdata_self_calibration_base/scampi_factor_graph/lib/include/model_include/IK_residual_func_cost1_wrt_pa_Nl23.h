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
 * Symbolic function: IK_residual_func_cost1_wrt_pa_Nl23
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
Eigen::Matrix<Scalar, 4, 3> IkResidualFuncCost1WrtPaNl23(
    const Scalar fh1, const Scalar fv1, const sym::Rot3<Scalar>& DeltaRot,
    const Eigen::Matrix<Scalar, 3, 1>& position_vector, const Eigen::Matrix<Scalar, 4, 1>& encoder,
    const Eigen::Matrix<Scalar, 3, 1>& p_a, const Eigen::Matrix<Scalar, 3, 1>& p_b,
    const Eigen::Matrix<Scalar, 3, 1>& p_c, const Eigen::Matrix<Scalar, 3, 1>& p_d,
    const sym::Rot3<Scalar>& Rot_init, const Scalar epsilon) {
  // Total ops: 1056

  // Unused inputs
  (void)encoder;
  (void)epsilon;

  // Input arrays
  const Eigen::Matrix<Scalar, 4, 1>& _DeltaRot = DeltaRot.Data();
  const Eigen::Matrix<Scalar, 4, 1>& _Rot_init = Rot_init.Data();

  // Intermediate terms (330)
  const Scalar _tmp0 = _DeltaRot[0] * _Rot_init[3] - _DeltaRot[1] * _Rot_init[2] +
                       _DeltaRot[2] * _Rot_init[1] + _DeltaRot[3] * _Rot_init[0];
  const Scalar _tmp1 = -2 * std::pow(_tmp0, Scalar(2));
  const Scalar _tmp2 = -_DeltaRot[0] * _Rot_init[1] + _DeltaRot[1] * _Rot_init[0] +
                       _DeltaRot[2] * _Rot_init[3] + _DeltaRot[3] * _Rot_init[2];
  const Scalar _tmp3 = 1 - 2 * std::pow(_tmp2, Scalar(2));
  const Scalar _tmp4 = Scalar(0.20999999999999999) * _tmp1 + Scalar(0.20999999999999999) * _tmp3;
  const Scalar _tmp5 = _DeltaRot[0] * _Rot_init[2] + _DeltaRot[1] * _Rot_init[3] -
                       _DeltaRot[2] * _Rot_init[0] + _DeltaRot[3] * _Rot_init[1];
  const Scalar _tmp6 = 2 * _tmp2;
  const Scalar _tmp7 = _tmp5 * _tmp6;
  const Scalar _tmp8 = -_DeltaRot[0] * _Rot_init[0] - _DeltaRot[1] * _Rot_init[1] -
                       _DeltaRot[2] * _Rot_init[2] + _DeltaRot[3] * _Rot_init[3];
  const Scalar _tmp9 = 2 * _tmp0;
  const Scalar _tmp10 = _tmp8 * _tmp9;
  const Scalar _tmp11 = -_tmp10 + _tmp7;
  const Scalar _tmp12 = -Scalar(0.010999999999999999) * _tmp11;
  const Scalar _tmp13 = _tmp5 * _tmp9;
  const Scalar _tmp14 = _tmp6 * _tmp8;
  const Scalar _tmp15 = Scalar(0.20999999999999999) * _tmp13 + Scalar(0.20999999999999999) * _tmp14;
  const Scalar _tmp16 = _tmp12 + _tmp15;
  const Scalar _tmp17 = _tmp16 + _tmp4;
  const Scalar _tmp18 = Scalar(1.0) * _tmp17;
  const Scalar _tmp19 = -_tmp18;
  const Scalar _tmp20 = -_tmp4;
  const Scalar _tmp21 = _tmp12 - _tmp15;
  const Scalar _tmp22 = _tmp20 + _tmp21;
  const Scalar _tmp23 = _tmp19 + _tmp22;
  const Scalar _tmp24 = _tmp16 + _tmp20;
  const Scalar _tmp25 = Scalar(1.0) / (_tmp19 + _tmp24);
  const Scalar _tmp26 = Scalar(0.20999999999999999) * _tmp13 - Scalar(0.20999999999999999) * _tmp14;
  const Scalar _tmp27 = -_tmp26;
  const Scalar _tmp28 = _tmp0 * _tmp6;
  const Scalar _tmp29 = 2 * _tmp5 * _tmp8;
  const Scalar _tmp30 = _tmp28 + _tmp29;
  const Scalar _tmp31 = -Scalar(0.010999999999999999) * _tmp30;
  const Scalar _tmp32 = -2 * std::pow(_tmp5, Scalar(2));
  const Scalar _tmp33 = Scalar(0.20999999999999999) * _tmp3 + Scalar(0.20999999999999999) * _tmp32;
  const Scalar _tmp34 = _tmp31 + _tmp33;
  const Scalar _tmp35 = _tmp27 + _tmp34;
  const Scalar _tmp36 = _tmp26 + _tmp34;
  const Scalar _tmp37 = Scalar(1.0) * _tmp36;
  const Scalar _tmp38 = _tmp25 * (-_tmp35 + _tmp37);
  const Scalar _tmp39 = _tmp23 * _tmp38;
  const Scalar _tmp40 = _tmp31 - _tmp33;
  const Scalar _tmp41 = _tmp27 + _tmp40;
  const Scalar _tmp42 = Scalar(1.0) / (_tmp37 - _tmp39 - _tmp41);
  const Scalar _tmp43 = _tmp18 * _tmp38 + _tmp37;
  const Scalar _tmp44 = _tmp42 * _tmp43;
  const Scalar _tmp45 = 0;
  const Scalar _tmp46 = _tmp41 + position_vector(0, 0);
  const Scalar _tmp47 = _tmp46 - p_a(0, 0);
  const Scalar _tmp48 = std::pow(_tmp47, Scalar(2));
  const Scalar _tmp49 = _tmp22 + position_vector(1, 0);
  const Scalar _tmp50 = _tmp49 - p_a(1, 0);
  const Scalar _tmp51 = std::pow(_tmp50, Scalar(2));
  const Scalar _tmp52 = _tmp48 + _tmp51;
  const Scalar _tmp53 = std::pow(_tmp52, Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp54 = Scalar(0.20999999999999999) * _tmp10 + Scalar(0.20999999999999999) * _tmp7;
  const Scalar _tmp55 = -_tmp54;
  const Scalar _tmp56 = -Scalar(0.010999999999999999) * _tmp1 -
                        Scalar(0.010999999999999999) * _tmp32 + Scalar(-0.010999999999999999);
  const Scalar _tmp57 = Scalar(0.20999999999999999) * _tmp28 - Scalar(0.20999999999999999) * _tmp29;
  const Scalar _tmp58 = _tmp56 - _tmp57;
  const Scalar _tmp59 = _tmp55 + _tmp58;
  const Scalar _tmp60 = _tmp53 * _tmp59;
  const Scalar _tmp61 = _tmp56 + _tmp57;
  const Scalar _tmp62 = _tmp54 + _tmp61;
  const Scalar _tmp63 = _tmp53 * _tmp62;
  const Scalar _tmp64 = _tmp47 * _tmp63;
  const Scalar _tmp65 = _tmp50 * _tmp53;
  const Scalar _tmp66 = _tmp36 + position_vector(0, 0);
  const Scalar _tmp67 = _tmp66 - p_c(0, 0);
  const Scalar _tmp68 = Scalar(1.0) / (_tmp67);
  const Scalar _tmp69 = _tmp17 + position_vector(1, 0);
  const Scalar _tmp70 = _tmp69 - p_c(1, 0);
  const Scalar _tmp71 = _tmp68 * _tmp70;
  const Scalar _tmp72 = _tmp53 * _tmp71;
  const Scalar _tmp73 = _tmp47 * _tmp72 - _tmp65;
  const Scalar _tmp74 = _tmp55 + _tmp61;
  const Scalar _tmp75 = _tmp35 + position_vector(0, 0);
  const Scalar _tmp76 = _tmp75 - p_b(0, 0);
  const Scalar _tmp77 = _tmp24 + position_vector(1, 0);
  const Scalar _tmp78 = _tmp77 - p_b(1, 0);
  const Scalar _tmp79 = std::pow(Scalar(std::pow(_tmp76, Scalar(2)) + std::pow(_tmp78, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp80 = _tmp78 * _tmp79;
  const Scalar _tmp81 = _tmp76 * _tmp79;
  const Scalar _tmp82 = _tmp62 * _tmp71;
  const Scalar _tmp83 = _tmp74 * _tmp80 - _tmp81 * _tmp82;
  const Scalar _tmp84 = Scalar(1.0) / (_tmp71 * _tmp81 - _tmp80);
  const Scalar _tmp85 = _tmp83 * _tmp84;
  const Scalar _tmp86 = _tmp59 * _tmp65 - _tmp64 * _tmp71 - _tmp73 * _tmp85;
  const Scalar _tmp87 = _tmp62 * _tmp81 - _tmp74 * _tmp81;
  const Scalar _tmp88 = _tmp84 * _tmp87;
  const Scalar _tmp89 = -_tmp38 * _tmp86 - _tmp47 * _tmp60 + _tmp64 - _tmp73 * _tmp88;
  const Scalar _tmp90 = Scalar(1.0) / (_tmp89);
  const Scalar _tmp91 = _tmp86 * _tmp90;
  const Scalar _tmp92 = _tmp19 - _tmp23 * _tmp44 - _tmp45 * _tmp91;
  const Scalar _tmp93 = Scalar(1.0) * _tmp25;
  const Scalar _tmp94 = Scalar(1.0) * _tmp42;
  const Scalar _tmp95 = Scalar(333.54000000000002) - fv1;
  const Scalar _tmp96 = _tmp23 * _tmp25;
  const Scalar _tmp97 = _tmp26 + _tmp40;
  const Scalar _tmp98 = _tmp97 - p_d(0, 0) + position_vector(0, 0);
  const Scalar _tmp99 = _tmp21 + _tmp4;
  const Scalar _tmp100 = _tmp99 - p_d(1, 0) + position_vector(1, 0);
  const Scalar _tmp101 =
      std::pow(Scalar(std::pow(_tmp100, Scalar(2)) + std::pow(_tmp98, Scalar(2))),
               Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp102 = _tmp101 * _tmp98;
  const Scalar _tmp103 = fh1 * (_tmp54 + _tmp58);
  const Scalar _tmp104 = _tmp102 * _tmp103 + Scalar(40.024799999999999) * _tmp30 + _tmp97 * fv1;
  const Scalar _tmp105 =
      std::sqrt(Scalar(std::pow(_tmp67, Scalar(2)) + std::pow(_tmp70, Scalar(2))));
  const Scalar _tmp106 = Scalar(1.0) / (_tmp105);
  const Scalar _tmp107 = _tmp105 * _tmp68;
  const Scalar _tmp108 = _tmp107 * (-_tmp106 * _tmp17 * _tmp67 + _tmp106 * _tmp36 * _tmp70);
  const Scalar _tmp109 = _tmp47 * _tmp53;
  const Scalar _tmp110 = _tmp108 * _tmp81 + _tmp24 * _tmp81 - _tmp35 * _tmp80;
  const Scalar _tmp111 = _tmp110 * _tmp84;
  const Scalar _tmp112 = _tmp22 * _tmp53;
  const Scalar _tmp113 = _tmp108 * _tmp109 - _tmp111 * _tmp73 + _tmp112 * _tmp47 - _tmp41 * _tmp65;
  const Scalar _tmp114 = Scalar(1.0) / (_tmp113);
  const Scalar _tmp115 = _tmp71 * _tmp85 + _tmp82;
  const Scalar _tmp116 = -_tmp115 * _tmp38 - _tmp62 + _tmp71 * _tmp88;
  const Scalar _tmp117 = _tmp113 * _tmp90;
  const Scalar _tmp118 = -_tmp108 + _tmp111 * _tmp71 - _tmp116 * _tmp117;
  const Scalar _tmp119 = _tmp114 * _tmp118;
  const Scalar _tmp120 = _tmp119 * _tmp89;
  const Scalar _tmp121 = _tmp116 + _tmp120;
  const Scalar _tmp122 = _tmp121 * _tmp90;
  const Scalar _tmp123 = _tmp23 * _tmp42;
  const Scalar _tmp124 = _tmp115 + _tmp120 * _tmp123 - _tmp122 * _tmp86;
  const Scalar _tmp125 = _tmp102 * fh1;
  const Scalar _tmp126 = Scalar(1.0) * _tmp125;
  const Scalar _tmp127 = Scalar(1.0) * _tmp84;
  const Scalar _tmp128 = _tmp127 * _tmp83;
  const Scalar _tmp129 = -_tmp127 * _tmp87 + _tmp128 * _tmp38;
  const Scalar _tmp130 = -_tmp110 * _tmp127 - _tmp117 * _tmp129;
  const Scalar _tmp131 = _tmp114 * _tmp89;
  const Scalar _tmp132 = _tmp130 * _tmp131;
  const Scalar _tmp133 = _tmp129 + _tmp132;
  const Scalar _tmp134 = _tmp133 * _tmp90;
  const Scalar _tmp135 = _tmp123 * _tmp132 - _tmp128 - _tmp134 * _tmp86;
  const Scalar _tmp136 = _tmp100 * _tmp101;
  const Scalar _tmp137 = _tmp136 * fh1;
  const Scalar _tmp138 = Scalar(1.0) * _tmp137;
  const Scalar _tmp139 = -_tmp103 * _tmp136 - Scalar(40.024799999999999) * _tmp11 - _tmp99 * fv1;
  const Scalar _tmp140 = _tmp39 * _tmp94 + Scalar(1.0);
  const Scalar _tmp141 = _tmp38 * _tmp94;
  const Scalar _tmp142 = Scalar(1.0) * _tmp114;
  const Scalar _tmp143 = _tmp131 * _tmp94;
  const Scalar _tmp144 = -_tmp142 * _tmp86 + _tmp143 * _tmp23;
  const Scalar _tmp145 = fh1 * (-_tmp102 * _tmp99 + _tmp136 * _tmp97);
  const Scalar _tmp146 = Scalar(1.0) * _tmp145;
  const Scalar _tmp147 = Scalar(1.0) * _tmp104 * (_tmp94 * _tmp96 - _tmp94) +
                         _tmp126 * (_tmp120 * _tmp94 - _tmp124 * _tmp93) +
                         _tmp138 * (_tmp132 * _tmp94 - _tmp135 * _tmp93) +
                         Scalar(1.0) * _tmp139 * (-_tmp140 * _tmp93 + _tmp141) +
                         _tmp146 * (_tmp143 - _tmp144 * _tmp93) +
                         Scalar(1.0) * _tmp95 * (-_tmp43 * _tmp94 - _tmp92 * _tmp93 + Scalar(1.0));
  const Scalar _tmp148 = _tmp73 * _tmp90;
  const Scalar _tmp149 = _tmp81 * _tmp84;
  const Scalar _tmp150 = _tmp149 * _tmp45;
  const Scalar _tmp151 = _tmp45 * _tmp90;
  const Scalar _tmp152 = _tmp151 * _tmp53;
  const Scalar _tmp153 = _tmp107 * _tmp95;
  const Scalar _tmp154 = -_tmp134 * _tmp73 + Scalar(1.0);
  const Scalar _tmp155 = _tmp134 * _tmp53;
  const Scalar _tmp156 = _tmp107 * _tmp137;
  const Scalar _tmp157 = _tmp142 * _tmp53;
  const Scalar _tmp158 = _tmp114 * _tmp127;
  const Scalar _tmp159 = _tmp158 * _tmp81;
  const Scalar _tmp160 = _tmp107 * _tmp145;
  const Scalar _tmp161 = -_tmp122 * _tmp73 - _tmp71;
  const Scalar _tmp162 = _tmp122 * _tmp53;
  const Scalar _tmp163 = _tmp107 * _tmp125;
  const Scalar _tmp164 = -_tmp153 * (-_tmp148 * _tmp150 + _tmp152 * _tmp47) -
                         _tmp156 * (_tmp149 * _tmp154 + _tmp155 * _tmp47) -
                         _tmp160 * (_tmp157 * _tmp47 - _tmp159 * _tmp73) -
                         _tmp163 * (_tmp149 * _tmp161 + _tmp162 * _tmp47 + Scalar(1.0));
  const Scalar _tmp165 = std::pow(_tmp164, Scalar(-2));
  const Scalar _tmp166 =
      std::pow(Scalar(std::pow(_tmp147, Scalar(2)) * _tmp165 + 1), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp167 = Scalar(1.0) / (_tmp164);
  const Scalar _tmp168 = std::pow(_tmp52, Scalar(Scalar(-3) / Scalar(2)));
  const Scalar _tmp169 = _tmp168 * _tmp48;
  const Scalar _tmp170 = _tmp168 * _tmp47 * _tmp50;
  const Scalar _tmp171 = _tmp169 * _tmp71 - _tmp170 - _tmp72;
  const Scalar _tmp172 = _tmp108 * _tmp169 - _tmp108 * _tmp53 - _tmp111 * _tmp171 - _tmp112 +
                         _tmp169 * _tmp22 - _tmp170 * _tmp41;
  const Scalar _tmp173 = _tmp172 * _tmp90;
  const Scalar _tmp174 = _tmp170 * _tmp59;
  const Scalar _tmp175 = -_tmp169 * _tmp82 - _tmp171 * _tmp85 + _tmp174 + _tmp63 * _tmp71;
  const Scalar _tmp176 =
      -_tmp169 * _tmp59 + _tmp169 * _tmp62 - _tmp171 * _tmp88 - _tmp175 * _tmp38 + _tmp60 - _tmp63;
  const Scalar _tmp177 = std::pow(_tmp89, Scalar(-2));
  const Scalar _tmp178 = _tmp176 * _tmp177;
  const Scalar _tmp179 = _tmp113 * _tmp178;
  const Scalar _tmp180 = _tmp131 * (-_tmp129 * _tmp173 + _tmp129 * _tmp179);
  const Scalar _tmp181 = _tmp114 * _tmp130;
  const Scalar _tmp182 = _tmp176 * _tmp181;
  const Scalar _tmp183 = std::pow(_tmp113, Scalar(-2));
  const Scalar _tmp184 = _tmp172 * _tmp183;
  const Scalar _tmp185 = _tmp130 * _tmp89;
  const Scalar _tmp186 = _tmp184 * _tmp185;
  const Scalar _tmp187 = _tmp180 + _tmp182 - _tmp186;
  const Scalar _tmp188 = _tmp133 * _tmp178;
  const Scalar _tmp189 = _tmp123 * _tmp180 + _tmp123 * _tmp182 - _tmp123 * _tmp186 -
                         _tmp134 * _tmp175 - _tmp187 * _tmp91 + _tmp188 * _tmp86;
  const Scalar _tmp190 = _tmp89 * _tmp94;
  const Scalar _tmp191 = _tmp184 * _tmp190;
  const Scalar _tmp192 = _tmp114 * _tmp94;
  const Scalar _tmp193 = _tmp176 * _tmp192;
  const Scalar _tmp194 = Scalar(1.0) * _tmp86;
  const Scalar _tmp195 =
      -_tmp142 * _tmp175 + _tmp184 * _tmp194 - _tmp191 * _tmp23 + _tmp193 * _tmp23;
  const Scalar _tmp196 = _tmp118 * _tmp89;
  const Scalar _tmp197 = _tmp184 * _tmp196;
  const Scalar _tmp198 = _tmp121 * _tmp178;
  const Scalar _tmp199 = _tmp119 * _tmp176;
  const Scalar _tmp200 = _tmp131 * (-_tmp116 * _tmp173 + _tmp116 * _tmp179);
  const Scalar _tmp201 = -_tmp197 + _tmp199 + _tmp200;
  const Scalar _tmp202 = -_tmp122 * _tmp175 - _tmp123 * _tmp197 + _tmp123 * _tmp199 +
                         _tmp123 * _tmp200 + _tmp198 * _tmp86 - _tmp201 * _tmp91;
  const Scalar _tmp203 = _tmp45 * _tmp86;
  const Scalar _tmp204 = -_tmp151 * _tmp175 + _tmp178 * _tmp203;
  const Scalar _tmp205 = _tmp93 * _tmp95;
  const Scalar _tmp206 = _tmp127 * _tmp73;
  const Scalar _tmp207 = _tmp206 * _tmp81;
  const Scalar _tmp208 = Scalar(1.0) * _tmp109;
  const Scalar _tmp209 = -_tmp122 * _tmp171 - _tmp148 * _tmp201 + _tmp198 * _tmp73;
  const Scalar _tmp210 = _tmp201 * _tmp90;
  const Scalar _tmp211 = _tmp149 * _tmp151;
  const Scalar _tmp212 = _tmp150 * _tmp73;
  const Scalar _tmp213 = _tmp109 * _tmp45;
  const Scalar _tmp214 = _tmp109 * _tmp90;
  const Scalar _tmp215 = -_tmp134 * _tmp171 - _tmp148 * _tmp187 + _tmp188 * _tmp73;
  const Scalar _tmp216 = -_tmp153 * (_tmp151 * _tmp169 - _tmp152 - _tmp171 * _tmp211 +
                                     _tmp178 * _tmp212 - _tmp178 * _tmp213) -
                         _tmp156 * (-_tmp109 * _tmp188 + _tmp134 * _tmp169 + _tmp149 * _tmp215 -
                                    _tmp155 + _tmp187 * _tmp214) -
                         _tmp160 * (_tmp142 * _tmp169 - _tmp157 - _tmp159 * _tmp171 +
                                    _tmp184 * _tmp207 - _tmp184 * _tmp208) -
                         _tmp163 * (-_tmp109 * _tmp198 + _tmp109 * _tmp210 + _tmp122 * _tmp169 +
                                    _tmp149 * _tmp209 - _tmp162);
  const Scalar _tmp217 = _tmp147 * _tmp165;
  const Scalar _tmp218 =
      _tmp166 *
      (_tmp167 *
           (_tmp126 * (-_tmp197 * _tmp94 + _tmp199 * _tmp94 + _tmp200 * _tmp94 - _tmp202 * _tmp93) +
            _tmp138 * (_tmp180 * _tmp94 + _tmp182 * _tmp94 - _tmp186 * _tmp94 - _tmp189 * _tmp93) +
            _tmp146 * (-_tmp191 + _tmp193 - _tmp195 * _tmp93) - _tmp204 * _tmp205) -
       _tmp216 * _tmp217);
  const Scalar _tmp219 = std::asinh(_tmp147 * _tmp167);
  const Scalar _tmp220 = Scalar(1.0) * _tmp219;
  const Scalar _tmp221 = Scalar(1.0) * std::sinh(_tmp220);
  const Scalar _tmp222 = Scalar(1.4083112389913199) * _tmp164;
  const Scalar _tmp223 =
      -_tmp219 * _tmp222 - std::sqrt(Scalar(std::pow(Scalar(-_tmp66 + p_c(0, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp69 + p_c(1, 0)), Scalar(2))));
  const Scalar _tmp224 = Scalar(0.71007031138673404) * _tmp165;
  const Scalar _tmp225 = _tmp223 * _tmp224;
  const Scalar _tmp226 = Scalar(1.4083112389913199) * _tmp219;
  const Scalar _tmp227 = Scalar(0.71007031138673404) * _tmp167;
  const Scalar _tmp228 = _tmp223 * _tmp227;
  const Scalar _tmp229 = std::sinh(_tmp228);
  const Scalar _tmp230 = _tmp224 * p_c(2, 0);
  const Scalar _tmp231 = Scalar(1.4083112389913199) * _tmp227 * p_c(2, 0) +
                         Scalar(1.4083112389913199) * std::cosh(_tmp220) -
                         Scalar(1.4083112389913199) * std::cosh(_tmp228);
  const Scalar _tmp232 = _tmp25 * _tmp95;
  const Scalar _tmp233 = _tmp137 * _tmp25;
  const Scalar _tmp234 = _tmp145 * _tmp25;
  const Scalar _tmp235 = _tmp125 * _tmp25;
  const Scalar _tmp236 = _tmp104 * _tmp94;
  const Scalar _tmp237 = _tmp124 * _tmp235 + _tmp135 * _tmp233 + _tmp139 * _tmp140 * _tmp25 +
                         _tmp144 * _tmp234 + _tmp232 * _tmp92 - _tmp236 * _tmp96;
  const Scalar _tmp238 = _tmp137 * _tmp84;
  const Scalar _tmp239 = _tmp125 * _tmp84;
  const Scalar _tmp240 = _tmp45 * _tmp95;
  const Scalar _tmp241 = _tmp145 * _tmp158;
  const Scalar _tmp242 =
      -_tmp148 * _tmp240 * _tmp84 + _tmp154 * _tmp238 + _tmp161 * _tmp239 - _tmp241 * _tmp73;
  const Scalar _tmp243 = Scalar(1.0) / (_tmp242);
  const Scalar _tmp244 = std::asinh(_tmp237 * _tmp243);
  const Scalar _tmp245 = Scalar(1.0) * _tmp244;
  const Scalar _tmp246 = Scalar(0.71007031138673404) * _tmp243;
  const Scalar _tmp247 = Scalar(1.4083112389913199) * _tmp242;
  const Scalar _tmp248 =
      -_tmp244 * _tmp247 - std::sqrt(Scalar(std::pow(Scalar(-_tmp75 + p_b(0, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp77 + p_b(1, 0)), Scalar(2))));
  const Scalar _tmp249 = _tmp246 * _tmp248;
  const Scalar _tmp250 = _tmp246 * p_b(2, 0) + std::cosh(_tmp245) - std::cosh(_tmp249);
  const Scalar _tmp251 = _tmp151 * _tmp95;
  const Scalar _tmp252 = _tmp251 * _tmp84;
  const Scalar _tmp253 = _tmp178 * _tmp240;
  const Scalar _tmp254 = _tmp73 * _tmp84;
  const Scalar _tmp255 = _tmp145 * _tmp206;
  const Scalar _tmp256 = -_tmp171 * _tmp241 - _tmp171 * _tmp252 + _tmp184 * _tmp255 +
                         _tmp209 * _tmp239 + _tmp215 * _tmp238 + _tmp253 * _tmp254;
  const Scalar _tmp257 = Scalar(1.4083112389913199) * _tmp256;
  const Scalar _tmp258 = std::pow(_tmp242, Scalar(-2));
  const Scalar _tmp259 = Scalar(0.71007031138673404) * _tmp258;
  const Scalar _tmp260 = _tmp256 * _tmp259;
  const Scalar _tmp261 =
      std::pow(Scalar(std::pow(_tmp237, Scalar(2)) * _tmp258 + 1), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp262 = _tmp237 * _tmp258;
  const Scalar _tmp263 =
      _tmp261 *
      (_tmp243 * (_tmp189 * _tmp233 + _tmp195 * _tmp234 + _tmp202 * _tmp235 + _tmp204 * _tmp232) -
       _tmp256 * _tmp262);
  const Scalar _tmp264 = std::sinh(_tmp249);
  const Scalar _tmp265 = Scalar(1.0) * std::sinh(_tmp245);
  const Scalar _tmp266 = _tmp137 * _tmp90;
  const Scalar _tmp267 = -_tmp125 * _tmp198 + _tmp125 * _tmp210 - _tmp137 * _tmp188 -
                         _tmp146 * _tmp184 + _tmp187 * _tmp266 - _tmp253;
  const Scalar _tmp268 = _tmp122 * _tmp125 + _tmp134 * _tmp137 + _tmp142 * _tmp145 + _tmp251;
  const Scalar _tmp269 = Scalar(1.0) / (_tmp268);
  const Scalar _tmp270 = Scalar(0.71007031138673404) * _tmp269;
  const Scalar _tmp271 = _tmp125 * _tmp42;
  const Scalar _tmp272 = _tmp137 * _tmp42;
  const Scalar _tmp273 = -_tmp120 * _tmp271 - _tmp132 * _tmp272 - _tmp139 * _tmp141 -
                         _tmp143 * _tmp145 + _tmp236 + _tmp44 * _tmp95;
  const Scalar _tmp274 = std::asinh(_tmp269 * _tmp273);
  const Scalar _tmp275 = Scalar(1.4083112389913199) * _tmp268;
  const Scalar _tmp276 = -_tmp49 + p_a(1, 0);
  const Scalar _tmp277 = -_tmp46 + p_a(0, 0);
  const Scalar _tmp278 =
      std::sqrt(Scalar(std::pow(_tmp276, Scalar(2)) + std::pow(_tmp277, Scalar(2))));
  const Scalar _tmp279 = -_tmp274 * _tmp275 - _tmp278;
  const Scalar _tmp280 = _tmp270 * _tmp279;
  const Scalar _tmp281 = Scalar(1.0) * _tmp274;
  const Scalar _tmp282 = Scalar(1.4083112389913199) * _tmp270 * p_a(2, 0) -
                         Scalar(1.4083112389913199) * std::cosh(_tmp280) +
                         Scalar(1.4083112389913199) * std::cosh(_tmp281);
  const Scalar _tmp283 = std::pow(_tmp268, Scalar(-2));
  const Scalar _tmp284 = Scalar(0.71007031138673404) * _tmp283;
  const Scalar _tmp285 = _tmp284 * p_a(2, 0);
  const Scalar _tmp286 = Scalar(1.0) / (_tmp278);
  const Scalar _tmp287 =
      std::pow(Scalar(std::pow(_tmp273, Scalar(2)) * _tmp283 + 1), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp288 = _tmp273 * _tmp283;
  const Scalar _tmp289 =
      _tmp287 *
      (-_tmp267 * _tmp288 +
       _tmp269 * (_tmp145 * _tmp191 - _tmp145 * _tmp193 - _tmp180 * _tmp272 - _tmp182 * _tmp272 +
                  _tmp186 * _tmp272 + _tmp197 * _tmp271 - _tmp199 * _tmp271 - _tmp200 * _tmp271));
  const Scalar _tmp290 = Scalar(1.4083112389913199) * _tmp274;
  const Scalar _tmp291 = _tmp279 * _tmp284;
  const Scalar _tmp292 = std::sinh(_tmp280);
  const Scalar _tmp293 = Scalar(1.0) * std::sinh(_tmp281);
  const Scalar _tmp294 = _tmp168 * _tmp51;
  const Scalar _tmp295 = _tmp170 * _tmp71 - _tmp294 + _tmp53;
  const Scalar _tmp296 =
      _tmp108 * _tmp170 - _tmp111 * _tmp295 + _tmp170 * _tmp22 - _tmp294 * _tmp41 + _tmp41 * _tmp53;
  const Scalar _tmp297 = _tmp183 * _tmp296;
  const Scalar _tmp298 = -_tmp170 * _tmp82 + _tmp294 * _tmp59 - _tmp295 * _tmp85 - _tmp60;
  const Scalar _tmp299 = _tmp170 * _tmp62 - _tmp174 - _tmp295 * _tmp88 - _tmp298 * _tmp38;
  const Scalar _tmp300 = _tmp177 * _tmp299;
  const Scalar _tmp301 = _tmp121 * _tmp300;
  const Scalar _tmp302 = _tmp113 * _tmp300;
  const Scalar _tmp303 = _tmp296 * _tmp90;
  const Scalar _tmp304 = _tmp131 * (_tmp116 * _tmp302 - _tmp116 * _tmp303);
  const Scalar _tmp305 = _tmp196 * _tmp297;
  const Scalar _tmp306 = _tmp119 * _tmp299;
  const Scalar _tmp307 = _tmp304 - _tmp305 + _tmp306;
  const Scalar _tmp308 = -_tmp122 * _tmp295 - _tmp148 * _tmp307 + _tmp301 * _tmp73;
  const Scalar _tmp309 = _tmp133 * _tmp300;
  const Scalar _tmp310 = _tmp181 * _tmp299;
  const Scalar _tmp311 = _tmp185 * _tmp297;
  const Scalar _tmp312 = _tmp131 * (_tmp129 * _tmp302 - _tmp129 * _tmp303);
  const Scalar _tmp313 = _tmp310 - _tmp311 + _tmp312;
  const Scalar _tmp314 = -_tmp134 * _tmp295 - _tmp148 * _tmp313 + _tmp309 * _tmp73;
  const Scalar _tmp315 =
      -_tmp153 * (_tmp151 * _tmp170 - _tmp211 * _tmp295 + _tmp212 * _tmp300 - _tmp213 * _tmp300) -
      _tmp156 * (-_tmp109 * _tmp309 + _tmp134 * _tmp170 + _tmp149 * _tmp314 + _tmp214 * _tmp313) -
      _tmp160 * (_tmp142 * _tmp170 - _tmp159 * _tmp295 + _tmp207 * _tmp297 - _tmp208 * _tmp297) -
      _tmp163 * (-_tmp109 * _tmp301 + _tmp122 * _tmp170 + _tmp149 * _tmp308 + _tmp214 * _tmp307);
  const Scalar _tmp316 = -_tmp122 * _tmp298 + _tmp123 * _tmp304 - _tmp123 * _tmp305 +
                         _tmp123 * _tmp306 + _tmp301 * _tmp86 - _tmp307 * _tmp91;
  const Scalar _tmp317 = -_tmp151 * _tmp298 + _tmp203 * _tmp300;
  const Scalar _tmp318 = _tmp192 * _tmp299;
  const Scalar _tmp319 = _tmp190 * _tmp297;
  const Scalar _tmp320 =
      -_tmp142 * _tmp298 + _tmp194 * _tmp297 + _tmp23 * _tmp318 - _tmp23 * _tmp319;
  const Scalar _tmp321 = _tmp123 * _tmp310 - _tmp123 * _tmp311 + _tmp123 * _tmp312 -
                         _tmp134 * _tmp298 + _tmp309 * _tmp86 - _tmp313 * _tmp91;
  const Scalar _tmp322 =
      _tmp166 *
      (_tmp167 *
           (_tmp126 * (_tmp304 * _tmp94 - _tmp305 * _tmp94 + _tmp306 * _tmp94 - _tmp316 * _tmp93) +
            _tmp138 * (_tmp310 * _tmp94 - _tmp311 * _tmp94 + _tmp312 * _tmp94 - _tmp321 * _tmp93) +
            _tmp146 * (_tmp318 - _tmp319 - _tmp320 * _tmp93) - _tmp205 * _tmp317) -
       _tmp217 * _tmp315);
  const Scalar _tmp323 = _tmp240 * _tmp300;
  const Scalar _tmp324 = _tmp238 * _tmp314 + _tmp239 * _tmp308 - _tmp241 * _tmp295 -
                         _tmp252 * _tmp295 + _tmp254 * _tmp323 + _tmp255 * _tmp297;
  const Scalar _tmp325 =
      _tmp261 *
      (_tmp243 * (_tmp232 * _tmp317 + _tmp233 * _tmp321 + _tmp234 * _tmp320 + _tmp235 * _tmp316) -
       _tmp262 * _tmp324);
  const Scalar _tmp326 = _tmp259 * _tmp324;
  const Scalar _tmp327 = Scalar(1.4083112389913199) * _tmp324;
  const Scalar _tmp328 = -_tmp125 * _tmp301 + _tmp125 * _tmp307 * _tmp90 - _tmp137 * _tmp309 -
                         _tmp146 * _tmp297 + _tmp266 * _tmp313 - _tmp323;
  const Scalar _tmp329 =
      _tmp287 *
      (_tmp269 * (-_tmp145 * _tmp318 + _tmp145 * _tmp319 - _tmp271 * _tmp304 + _tmp271 * _tmp305 -
                  _tmp271 * _tmp306 - _tmp272 * _tmp310 + _tmp272 * _tmp311 - _tmp272 * _tmp312) -
       _tmp288 * _tmp328);

  // Output terms (1)
  Eigen::Matrix<Scalar, 4, 3> _res;

  _res(0, 0) = 0;
  _res(1, 0) =
      -_tmp216 * _tmp231 -
      _tmp222 *
          (-_tmp216 * _tmp230 + _tmp218 * _tmp221 -
           _tmp229 * (-_tmp216 * _tmp225 + _tmp227 * (-_tmp216 * _tmp226 - _tmp218 * _tmp222)));
  _res(2, 0) =
      -_tmp247 *
          (-_tmp260 * p_b(2, 0) + _tmp263 * _tmp265 -
           _tmp264 * (_tmp246 * (-_tmp244 * _tmp257 - _tmp247 * _tmp263) - _tmp248 * _tmp260)) -
      _tmp250 * _tmp257;
  _res(3, 0) =
      -_tmp267 * _tmp282 -
      _tmp275 * (-_tmp267 * _tmp285 + _tmp289 * _tmp293 -
                 _tmp292 * (-_tmp267 * _tmp291 + _tmp270 * (-_tmp267 * _tmp290 - _tmp275 * _tmp289 -
                                                            _tmp277 * _tmp286)));
  _res(0, 1) = 0;
  _res(1, 1) =
      -_tmp222 *
          (_tmp221 * _tmp322 -
           _tmp229 * (-_tmp225 * _tmp315 + _tmp227 * (-_tmp222 * _tmp322 - _tmp226 * _tmp315)) -
           _tmp230 * _tmp315) -
      _tmp231 * _tmp315;
  _res(2, 1) =
      -_tmp247 *
          (-_tmp264 * (_tmp246 * (-_tmp244 * _tmp327 - _tmp247 * _tmp325) - _tmp248 * _tmp326) +
           _tmp265 * _tmp325 - _tmp326 * p_b(2, 0)) -
      _tmp250 * _tmp327;
  _res(3, 1) =
      -_tmp275 *
          (-_tmp285 * _tmp328 -
           _tmp292 * (_tmp270 * (-_tmp275 * _tmp329 - _tmp276 * _tmp286 - _tmp290 * _tmp328) -
                      _tmp291 * _tmp328) +
           _tmp293 * _tmp329) -
      _tmp282 * _tmp328;
  _res(0, 2) = 0;
  _res(1, 2) = 0;
  _res(2, 2) = 0;
  _res(3, 2) = Scalar(-1.0);

  return _res;
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym

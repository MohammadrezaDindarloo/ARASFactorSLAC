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
 * Symbolic function: IK_residual_func_cost1_wrt_pd_Nl20
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
Eigen::Matrix<Scalar, 4, 3> IkResidualFuncCost1WrtPdNl20(
    const Scalar fh1, const Scalar fv1, const sym::Rot3<Scalar>& DeltaRot,
    const Eigen::Matrix<Scalar, 3, 1>& position_vector, const Eigen::Matrix<Scalar, 4, 1>& encoder,
    const Eigen::Matrix<Scalar, 3, 1>& p_a, const Eigen::Matrix<Scalar, 3, 1>& p_b,
    const Eigen::Matrix<Scalar, 3, 1>& p_c, const Eigen::Matrix<Scalar, 3, 1>& p_d,
    const sym::Rot3<Scalar>& Rot_init, const Scalar epsilon) {
  // Total ops: 777

  // Unused inputs
  (void)encoder;
  (void)epsilon;

  // Input arrays
  const Eigen::Matrix<Scalar, 4, 1>& _DeltaRot = DeltaRot.Data();
  const Eigen::Matrix<Scalar, 4, 1>& _Rot_init = Rot_init.Data();

  // Intermediate terms (270)
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
  const Scalar _tmp6 = Scalar(0.20999999999999999) * _tmp2 - Scalar(0.20999999999999999) * _tmp5;
  const Scalar _tmp7 = 2 * _tmp3;
  const Scalar _tmp8 = _tmp1 * _tmp7;
  const Scalar _tmp9 = _tmp0 * _tmp4;
  const Scalar _tmp10 = _tmp8 + _tmp9;
  const Scalar _tmp11 = -Scalar(0.010999999999999999) * _tmp10;
  const Scalar _tmp12 = -2 * std::pow(_tmp0, Scalar(2));
  const Scalar _tmp13 = 1 - 2 * std::pow(_tmp3, Scalar(2));
  const Scalar _tmp14 = Scalar(0.20999999999999999) * _tmp12 + Scalar(0.20999999999999999) * _tmp13;
  const Scalar _tmp15 = _tmp11 - _tmp14;
  const Scalar _tmp16 = _tmp15 + _tmp6;
  const Scalar _tmp17 = _tmp16 + position_vector(0, 0);
  const Scalar _tmp18 = -_tmp17 + p_d(0, 0);
  const Scalar _tmp19 = Scalar(0.20999999999999999) * _tmp2 + Scalar(0.20999999999999999) * _tmp5;
  const Scalar _tmp20 = -_tmp19;
  const Scalar _tmp21 = -2 * std::pow(_tmp1, Scalar(2));
  const Scalar _tmp22 = Scalar(0.20999999999999999) * _tmp13 + Scalar(0.20999999999999999) * _tmp21;
  const Scalar _tmp23 = _tmp0 * _tmp7;
  const Scalar _tmp24 = _tmp1 * _tmp4;
  const Scalar _tmp25 = _tmp23 - _tmp24;
  const Scalar _tmp26 = -Scalar(0.010999999999999999) * _tmp25;
  const Scalar _tmp27 = _tmp22 + _tmp26;
  const Scalar _tmp28 = _tmp20 + _tmp27;
  const Scalar _tmp29 = _tmp28 + position_vector(1, 0);
  const Scalar _tmp30 = -_tmp29 + p_d(1, 0);
  const Scalar _tmp31 =
      std::sqrt(Scalar(std::pow(_tmp18, Scalar(2)) + std::pow(_tmp30, Scalar(2))));
  const Scalar _tmp32 = Scalar(1.0) / (fh1);
  const Scalar _tmp33 =
      Scalar(1.0) *
      std::sinh(Scalar(0.71007031138673404) * _tmp32 *
                (-_tmp31 - Scalar(1.4083112389913199) * fh1 * std::asinh(_tmp32 * fv1))) /
      _tmp31;
  const Scalar _tmp34 = _tmp29 - p_d(1, 0);
  const Scalar _tmp35 = _tmp11 + _tmp14;
  const Scalar _tmp36 = _tmp35 + _tmp6;
  const Scalar _tmp37 = _tmp36 + position_vector(0, 0);
  const Scalar _tmp38 = _tmp37 - p_c(0, 0);
  const Scalar _tmp39 = _tmp19 + _tmp27;
  const Scalar _tmp40 = _tmp39 + position_vector(1, 0);
  const Scalar _tmp41 = _tmp40 - p_c(1, 0);
  const Scalar _tmp42 = std::pow(Scalar(std::pow(_tmp38, Scalar(2)) + std::pow(_tmp41, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp43 = _tmp38 * _tmp42;
  const Scalar _tmp44 = -_tmp6;
  const Scalar _tmp45 = _tmp15 + _tmp44;
  const Scalar _tmp46 = _tmp45 + position_vector(0, 0);
  const Scalar _tmp47 = _tmp46 - p_a(0, 0);
  const Scalar _tmp48 = -_tmp22 + _tmp26;
  const Scalar _tmp49 = _tmp20 + _tmp48;
  const Scalar _tmp50 = _tmp49 + position_vector(1, 0);
  const Scalar _tmp51 = _tmp50 - p_a(1, 0);
  const Scalar _tmp52 = std::pow(Scalar(std::pow(_tmp47, Scalar(2)) + std::pow(_tmp51, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp53 = _tmp47 * _tmp52;
  const Scalar _tmp54 = Scalar(0.20999999999999999) * _tmp8 - Scalar(0.20999999999999999) * _tmp9;
  const Scalar _tmp55 = -Scalar(0.010999999999999999) * _tmp12 -
                        Scalar(0.010999999999999999) * _tmp21 + Scalar(-0.010999999999999999);
  const Scalar _tmp56 = Scalar(0.20999999999999999) * _tmp23 + Scalar(0.20999999999999999) * _tmp24;
  const Scalar _tmp57 = _tmp55 - _tmp56;
  const Scalar _tmp58 = _tmp54 + _tmp57;
  const Scalar _tmp59 = _tmp35 + _tmp44;
  const Scalar _tmp60 = _tmp59 + position_vector(0, 0);
  const Scalar _tmp61 = _tmp60 - p_b(0, 0);
  const Scalar _tmp62 = Scalar(1.0) / (_tmp61);
  const Scalar _tmp63 = _tmp19 + _tmp48;
  const Scalar _tmp64 = _tmp63 + position_vector(1, 0);
  const Scalar _tmp65 = _tmp64 - p_b(1, 0);
  const Scalar _tmp66 = _tmp62 * _tmp65;
  const Scalar _tmp67 = _tmp58 * _tmp66;
  const Scalar _tmp68 = -_tmp54;
  const Scalar _tmp69 = _tmp57 + _tmp68;
  const Scalar _tmp70 = _tmp51 * _tmp52;
  const Scalar _tmp71 = -_tmp53 * _tmp67 + _tmp69 * _tmp70;
  const Scalar _tmp72 = Scalar(1.0) / (_tmp53 * _tmp66 - _tmp70);
  const Scalar _tmp73 = _tmp41 * _tmp42;
  const Scalar _tmp74 = _tmp43 * _tmp66 - _tmp73;
  const Scalar _tmp75 = _tmp72 * _tmp74;
  const Scalar _tmp76 = _tmp55 + _tmp56;
  const Scalar _tmp77 = _tmp54 + _tmp76;
  const Scalar _tmp78 = -_tmp43 * _tmp67 - _tmp71 * _tmp75 + _tmp73 * _tmp77;
  const Scalar _tmp79 = Scalar(1.0) * _tmp63;
  const Scalar _tmp80 = -_tmp79;
  const Scalar _tmp81 = Scalar(1.0) / (_tmp49 + _tmp80);
  const Scalar _tmp82 = Scalar(1.0) * _tmp59;
  const Scalar _tmp83 = _tmp81 * (-_tmp45 + _tmp82);
  const Scalar _tmp84 = _tmp53 * _tmp58 - _tmp53 * _tmp69;
  const Scalar _tmp85 = _tmp43 * _tmp58 - _tmp43 * _tmp77 - _tmp75 * _tmp84 - _tmp78 * _tmp83;
  const Scalar _tmp86 = Scalar(1.0) / (_tmp85);
  const Scalar _tmp87 = Scalar(1.0) * _tmp72;
  const Scalar _tmp88 = _tmp71 * _tmp87;
  const Scalar _tmp89 = _tmp83 * _tmp88 - _tmp84 * _tmp87;
  const Scalar _tmp90 =
      std::sqrt(Scalar(std::pow(_tmp61, Scalar(2)) + std::pow(_tmp65, Scalar(2))));
  const Scalar _tmp91 = Scalar(1.0) / (_tmp90);
  const Scalar _tmp92 = _tmp62 * _tmp90;
  const Scalar _tmp93 = _tmp92 * (_tmp59 * _tmp65 * _tmp91 - _tmp61 * _tmp63 * _tmp91);
  const Scalar _tmp94 = -_tmp45 * _tmp70 + _tmp49 * _tmp53 + _tmp53 * _tmp93;
  const Scalar _tmp95 = -_tmp36 * _tmp73 + _tmp39 * _tmp43 + _tmp43 * _tmp93 - _tmp75 * _tmp94;
  const Scalar _tmp96 = _tmp86 * _tmp95;
  const Scalar _tmp97 = Scalar(1.0) / (_tmp95);
  const Scalar _tmp98 = _tmp85 * _tmp97;
  const Scalar _tmp99 = _tmp98 * (-_tmp87 * _tmp94 - _tmp89 * _tmp96);
  const Scalar _tmp100 = _tmp89 + _tmp99;
  const Scalar _tmp101 = _tmp100 * _tmp86;
  const Scalar _tmp102 = _tmp74 * _tmp86;
  const Scalar _tmp103 = _tmp72 * (-_tmp100 * _tmp102 + Scalar(1.0));
  const Scalar _tmp104 = _tmp101 * _tmp43 + _tmp103 * _tmp53;
  const Scalar _tmp105 = _tmp17 - p_d(0, 0);
  const Scalar _tmp106 = std::pow(_tmp105, Scalar(2));
  const Scalar _tmp107 = std::pow(_tmp34, Scalar(2));
  const Scalar _tmp108 = _tmp106 + _tmp107;
  const Scalar _tmp109 = std::pow(_tmp108, Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp110 = _tmp109 * fh1;
  const Scalar _tmp111 = _tmp110 * _tmp92;
  const Scalar _tmp112 = _tmp104 * _tmp111;
  const Scalar _tmp113 = _tmp66 * _tmp72;
  const Scalar _tmp114 = _tmp113 * _tmp71 + _tmp67;
  const Scalar _tmp115 = _tmp113 * _tmp84 - _tmp114 * _tmp83 - _tmp58;
  const Scalar _tmp116 = _tmp98 * (_tmp113 * _tmp94 - _tmp115 * _tmp96 - _tmp93);
  const Scalar _tmp117 = _tmp115 + _tmp116;
  const Scalar _tmp118 = _tmp117 * _tmp86;
  const Scalar _tmp119 = _tmp72 * (-_tmp102 * _tmp117 - _tmp66);
  const Scalar _tmp120 = _tmp118 * _tmp43 + _tmp119 * _tmp53 + Scalar(1.0);
  const Scalar _tmp121 = _tmp111 * _tmp120;
  const Scalar _tmp122 = _tmp109 * _tmp16;
  const Scalar _tmp123 = _tmp109 * _tmp28;
  const Scalar _tmp124 = fh1 * (-_tmp105 * _tmp123 + _tmp122 * _tmp34);
  const Scalar _tmp125 = Scalar(1.0) * _tmp97;
  const Scalar _tmp126 = _tmp74 * _tmp87 * _tmp97;
  const Scalar _tmp127 = _tmp92 * (_tmp125 * _tmp43 - _tmp126 * _tmp53);
  const Scalar _tmp128 = Scalar(333.54000000000002) - fv1;
  const Scalar _tmp129 = _tmp79 * _tmp83 + _tmp82;
  const Scalar _tmp130 = 0;
  const Scalar _tmp131 = _tmp130 * _tmp86;
  const Scalar _tmp132 = -_tmp105 * _tmp121 - _tmp112 * _tmp34 - _tmp124 * _tmp127 -
                         _tmp128 * _tmp92 * (_tmp131 * _tmp43 - _tmp131 * _tmp53 * _tmp75);
  const Scalar _tmp133 = Scalar(1.0) / (_tmp132);
  const Scalar _tmp134 = _tmp39 + _tmp80;
  const Scalar _tmp135 = _tmp134 * _tmp83;
  const Scalar _tmp136 = Scalar(1.0) / (-_tmp135 - _tmp36 + _tmp82);
  const Scalar _tmp137 = Scalar(1.0) * _tmp136;
  const Scalar _tmp138 = _tmp78 * _tmp86;
  const Scalar _tmp139 = _tmp134 * _tmp136;
  const Scalar _tmp140 = -_tmp100 * _tmp138 + _tmp139 * _tmp99 - _tmp88;
  const Scalar _tmp141 = Scalar(1.0) * _tmp81;
  const Scalar _tmp142 = _tmp137 * _tmp99 - _tmp140 * _tmp141;
  const Scalar _tmp143 = Scalar(1.0) * _tmp110;
  const Scalar _tmp144 = _tmp142 * _tmp143;
  const Scalar _tmp145 = _tmp68 + _tmp76;
  const Scalar _tmp146 = _tmp110 * _tmp145;
  const Scalar _tmp147 = Scalar(40.024799999999999) * _tmp10 + _tmp105 * _tmp146 + _tmp16 * fv1;
  const Scalar _tmp148 = _tmp134 * _tmp81;
  const Scalar _tmp149 = _tmp137 * _tmp148;
  const Scalar _tmp150 = -Scalar(1.0) * _tmp137 + Scalar(1.0) * _tmp149;
  const Scalar _tmp151 = _tmp114 + _tmp116 * _tmp139 - _tmp117 * _tmp138;
  const Scalar _tmp152 = _tmp116 * _tmp137 - _tmp141 * _tmp151;
  const Scalar _tmp153 = _tmp143 * _tmp152;
  const Scalar _tmp154 = -_tmp146 * _tmp34 - Scalar(40.024799999999999) * _tmp25 - _tmp28 * fv1;
  const Scalar _tmp155 = _tmp137 * _tmp83;
  const Scalar _tmp156 = _tmp81 * (_tmp135 * _tmp137 + Scalar(1.0));
  const Scalar _tmp157 = Scalar(1.0) * _tmp155 - Scalar(1.0) * _tmp156;
  const Scalar _tmp158 = _tmp129 * _tmp136;
  const Scalar _tmp159 = -_tmp130 * _tmp138 - _tmp134 * _tmp158 + _tmp80;
  const Scalar _tmp160 = _tmp137 * _tmp98;
  const Scalar _tmp161 = -_tmp125 * _tmp78 + _tmp134 * _tmp160;
  const Scalar _tmp162 = -Scalar(1.0) * _tmp141 * _tmp161 + Scalar(1.0) * _tmp160;
  const Scalar _tmp163 =
      _tmp105 * _tmp153 + _tmp124 * _tmp162 +
      Scalar(1.0) * _tmp128 * (-_tmp129 * _tmp137 - _tmp141 * _tmp159 + Scalar(1.0)) +
      _tmp144 * _tmp34 + _tmp147 * _tmp150 + _tmp154 * _tmp157;
  const Scalar _tmp164 = std::asinh(_tmp133 * _tmp163);
  const Scalar _tmp165 = Scalar(1.0) * _tmp164;
  const Scalar _tmp166 = Scalar(1.0) * std::sinh(_tmp165);
  const Scalar _tmp167 = std::pow(_tmp132, Scalar(-2));
  const Scalar _tmp168 =
      std::pow(Scalar(std::pow(_tmp163, Scalar(2)) * _tmp167 + 1), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp169 = _tmp145 * fh1;
  const Scalar _tmp170 = std::pow(_tmp108, Scalar(Scalar(-3) / Scalar(2)));
  const Scalar _tmp171 = _tmp105 * _tmp34;
  const Scalar _tmp172 = _tmp170 * _tmp171;
  const Scalar _tmp173 = _tmp169 * _tmp172;
  const Scalar _tmp174 = _tmp106 * _tmp170;
  const Scalar _tmp175 = _tmp174 * fh1;
  const Scalar _tmp176 = Scalar(1.0) * _tmp152;
  const Scalar _tmp177 = _tmp170 * _tmp28;
  const Scalar _tmp178 = fh1 * (-_tmp106 * _tmp177 + _tmp123 + _tmp16 * _tmp172);
  const Scalar _tmp179 = Scalar(1.0) * _tmp142;
  const Scalar _tmp180 = _tmp172 * fh1;
  const Scalar _tmp181 = -_tmp146 + _tmp169 * _tmp174;
  const Scalar _tmp182 = _tmp104 * _tmp92;
  const Scalar _tmp183 = _tmp120 * _tmp92;
  const Scalar _tmp184 = _tmp121 - _tmp127 * _tmp178 - _tmp175 * _tmp183 - _tmp180 * _tmp182;
  const Scalar _tmp185 = _tmp163 * _tmp167;
  const Scalar _tmp186 =
      _tmp168 * (_tmp133 * (_tmp150 * _tmp181 - _tmp153 - _tmp157 * _tmp173 + _tmp162 * _tmp178 +
                            _tmp175 * _tmp176 + _tmp179 * _tmp180) -
                 _tmp184 * _tmp185);
  const Scalar _tmp187 = Scalar(1.4083112389913199) * _tmp132;
  const Scalar _tmp188 = Scalar(1.4083112389913199) * _tmp164;
  const Scalar _tmp189 = Scalar(0.71007031138673404) * _tmp133;
  const Scalar _tmp190 =
      -_tmp164 * _tmp187 - std::sqrt(Scalar(std::pow(Scalar(-_tmp60 + p_b(0, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp64 + p_b(1, 0)), Scalar(2))));
  const Scalar _tmp191 = Scalar(0.71007031138673404) * _tmp167;
  const Scalar _tmp192 = _tmp190 * _tmp191;
  const Scalar _tmp193 = _tmp189 * _tmp190;
  const Scalar _tmp194 = std::sinh(_tmp193);
  const Scalar _tmp195 = _tmp191 * p_b(2, 0);
  const Scalar _tmp196 = Scalar(1.4083112389913199) * _tmp189 * p_b(2, 0) +
                         Scalar(1.4083112389913199) * std::cosh(_tmp165) -
                         Scalar(1.4083112389913199) * std::cosh(_tmp193);
  const Scalar _tmp197 = _tmp161 * _tmp81;
  const Scalar _tmp198 = _tmp137 * _tmp147;
  const Scalar _tmp199 = _tmp110 * _tmp81;
  const Scalar _tmp200 = _tmp151 * _tmp199;
  const Scalar _tmp201 = _tmp140 * _tmp199;
  const Scalar _tmp202 = _tmp105 * _tmp200 + _tmp124 * _tmp197 + _tmp128 * _tmp159 * _tmp81 -
                         _tmp148 * _tmp198 + _tmp154 * _tmp156 + _tmp201 * _tmp34;
  const Scalar _tmp203 = _tmp128 * _tmp131;
  const Scalar _tmp204 = _tmp110 * _tmp119;
  const Scalar _tmp205 = _tmp103 * _tmp110;
  const Scalar _tmp206 =
      _tmp105 * _tmp204 - _tmp124 * _tmp126 - _tmp203 * _tmp75 + _tmp205 * _tmp34;
  const Scalar _tmp207 = Scalar(1.0) / (_tmp206);
  const Scalar _tmp208 = std::asinh(_tmp202 * _tmp207);
  const Scalar _tmp209 = Scalar(1.0) * _tmp208;
  const Scalar _tmp210 = Scalar(1.0) * std::sinh(_tmp209);
  const Scalar _tmp211 = std::pow(_tmp206, Scalar(-2));
  const Scalar _tmp212 =
      std::pow(Scalar(std::pow(_tmp202, Scalar(2)) * _tmp211 + 1), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp213 = _tmp103 * _tmp180 + _tmp119 * _tmp175 - _tmp126 * _tmp178 - _tmp204;
  const Scalar _tmp214 = _tmp202 * _tmp211;
  const Scalar _tmp215 = _tmp81 * fh1;
  const Scalar _tmp216 = _tmp140 * _tmp215;
  const Scalar _tmp217 = _tmp137 * _tmp181;
  const Scalar _tmp218 = _tmp151 * _tmp215;
  const Scalar _tmp219 =
      _tmp212 * (_tmp207 * (-_tmp148 * _tmp217 - _tmp156 * _tmp173 + _tmp172 * _tmp216 +
                            _tmp174 * _tmp218 + _tmp178 * _tmp197 - _tmp200) -
                 _tmp213 * _tmp214);
  const Scalar _tmp220 = Scalar(1.4083112389913199) * _tmp208;
  const Scalar _tmp221 =
      -_tmp206 * _tmp220 - std::sqrt(Scalar(std::pow(Scalar(-_tmp46 + p_a(0, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp50 + p_a(1, 0)), Scalar(2))));
  const Scalar _tmp222 = Scalar(0.71007031138673404) * _tmp207;
  const Scalar _tmp223 = _tmp221 * _tmp222;
  const Scalar _tmp224 = std::sinh(_tmp223);
  const Scalar _tmp225 = Scalar(1.4083112389913199) * _tmp206;
  const Scalar _tmp226 = Scalar(1.4083112389913199) * _tmp213;
  const Scalar _tmp227 = Scalar(0.71007031138673404) * _tmp211;
  const Scalar _tmp228 = _tmp213 * _tmp227;
  const Scalar _tmp229 = _tmp222 * p_a(2, 0) + std::cosh(_tmp209) - std::cosh(_tmp223);
  const Scalar _tmp230 = _tmp110 * _tmp118;
  const Scalar _tmp231 = _tmp101 * _tmp180 + _tmp118 * _tmp175 + _tmp125 * _tmp178 - _tmp230;
  const Scalar _tmp232 = _tmp110 * _tmp136;
  const Scalar _tmp233 = _tmp232 * _tmp99;
  const Scalar _tmp234 = _tmp116 * _tmp232;
  const Scalar _tmp235 = -_tmp105 * _tmp234 - _tmp124 * _tmp160 + _tmp128 * _tmp158 -
                         _tmp154 * _tmp155 + _tmp198 - _tmp233 * _tmp34;
  const Scalar _tmp236 = _tmp101 * _tmp110;
  const Scalar _tmp237 = _tmp105 * _tmp230 + _tmp124 * _tmp125 + _tmp203 + _tmp236 * _tmp34;
  const Scalar _tmp238 = Scalar(1.0) / (_tmp237);
  const Scalar _tmp239 = std::asinh(_tmp235 * _tmp238);
  const Scalar _tmp240 = Scalar(1.0) * _tmp239;
  const Scalar _tmp241 = Scalar(0.71007031138673404) * _tmp238;
  const Scalar _tmp242 = Scalar(1.4083112389913199) * _tmp239;
  const Scalar _tmp243 =
      -_tmp237 * _tmp242 - std::sqrt(Scalar(std::pow(Scalar(-_tmp37 + p_c(0, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp40 + p_c(1, 0)), Scalar(2))));
  const Scalar _tmp244 = _tmp241 * _tmp243;
  const Scalar _tmp245 = Scalar(1.4083112389913199) * _tmp241 * p_c(2, 0) +
                         Scalar(1.4083112389913199) * std::cosh(_tmp240) -
                         Scalar(1.4083112389913199) * std::cosh(_tmp244);
  const Scalar _tmp246 = std::pow(_tmp237, Scalar(-2));
  const Scalar _tmp247 = Scalar(0.71007031138673404) * _tmp246;
  const Scalar _tmp248 = _tmp231 * _tmp247;
  const Scalar _tmp249 = _tmp235 * _tmp246;
  const Scalar _tmp250 = _tmp136 * _tmp99;
  const Scalar _tmp251 = _tmp116 * _tmp136;
  const Scalar _tmp252 =
      -_tmp231 * _tmp249 + _tmp238 * (_tmp155 * _tmp173 - _tmp160 * _tmp178 - _tmp175 * _tmp251 -
                                      _tmp180 * _tmp250 + _tmp217 + _tmp234);
  const Scalar _tmp253 =
      std::pow(Scalar(std::pow(_tmp235, Scalar(2)) * _tmp246 + 1), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp254 = Scalar(1.4083112389913199) * _tmp237;
  const Scalar _tmp255 = _tmp253 * _tmp254;
  const Scalar _tmp256 = std::sinh(_tmp244);
  const Scalar _tmp257 = Scalar(1.0) * _tmp253 * std::sinh(_tmp240);
  const Scalar _tmp258 = _tmp107 * _tmp170;
  const Scalar _tmp259 = _tmp258 * fh1;
  const Scalar _tmp260 = fh1 * (-_tmp122 + _tmp16 * _tmp258 - _tmp171 * _tmp177);
  const Scalar _tmp261 = _tmp112 - _tmp127 * _tmp260 - _tmp180 * _tmp183 - _tmp182 * _tmp259;
  const Scalar _tmp262 = _tmp146 - _tmp169 * _tmp258;
  const Scalar _tmp263 =
      _tmp168 * (_tmp133 * (-_tmp144 + _tmp150 * _tmp173 + _tmp157 * _tmp262 + _tmp162 * _tmp260 +
                            _tmp176 * _tmp180 + _tmp179 * _tmp259) -
                 _tmp185 * _tmp261);
  const Scalar _tmp264 = _tmp103 * _tmp259 + _tmp119 * _tmp180 - _tmp126 * _tmp260 - _tmp205;
  const Scalar _tmp265 = _tmp227 * _tmp264;
  const Scalar _tmp266 =
      _tmp212 * (_tmp207 * (-_tmp149 * _tmp173 + _tmp156 * _tmp262 + _tmp172 * _tmp218 +
                            _tmp197 * _tmp260 - _tmp201 + _tmp216 * _tmp258) -
                 _tmp214 * _tmp264);
  const Scalar _tmp267 = _tmp101 * _tmp259 + _tmp118 * _tmp180 + _tmp125 * _tmp260 - _tmp236;
  const Scalar _tmp268 = _tmp247 * _tmp267;
  const Scalar _tmp269 = _tmp238 * (_tmp137 * _tmp173 - _tmp155 * _tmp262 - _tmp160 * _tmp260 -
                                    _tmp180 * _tmp251 + _tmp233 - _tmp250 * _tmp259) -
                         _tmp249 * _tmp267;

  // Output terms (1)
  Eigen::Matrix<Scalar, 4, 3> _res;

  _res(0, 0) = -_tmp18 * _tmp33;
  _res(1, 0) =
      -_tmp184 * _tmp196 -
      _tmp187 *
          (_tmp166 * _tmp186 - _tmp184 * _tmp195 -
           _tmp194 * (-_tmp184 * _tmp192 + _tmp189 * (-_tmp184 * _tmp188 - _tmp186 * _tmp187)));
  _res(2, 0) =
      -_tmp225 *
          (_tmp210 * _tmp219 -
           _tmp224 * (-_tmp221 * _tmp228 + _tmp222 * (-_tmp208 * _tmp226 - _tmp219 * _tmp225)) -
           _tmp228 * p_a(2, 0)) -
      _tmp226 * _tmp229;
  _res(3, 0) =
      -_tmp231 * _tmp245 -
      _tmp254 *
          (-_tmp248 * p_c(2, 0) + _tmp252 * _tmp257 -
           _tmp256 * (_tmp241 * (-_tmp231 * _tmp242 - _tmp252 * _tmp255) - _tmp243 * _tmp248));
  _res(0, 1) = -_tmp30 * _tmp33;
  _res(1, 1) =
      -_tmp187 *
          (_tmp166 * _tmp263 -
           _tmp194 * (_tmp189 * (-_tmp187 * _tmp263 - _tmp188 * _tmp261) - _tmp192 * _tmp261) -
           _tmp195 * _tmp261) -
      _tmp196 * _tmp261;
  _res(2, 1) =
      -_tmp225 *
          (_tmp210 * _tmp266 -
           _tmp224 * (-_tmp221 * _tmp265 + _tmp222 * (-_tmp220 * _tmp264 - _tmp225 * _tmp266)) -
           _tmp265 * p_a(2, 0)) -
      Scalar(1.4083112389913199) * _tmp229 * _tmp264;
  _res(3, 1) =
      -_tmp245 * _tmp267 -
      _tmp254 *
          (-_tmp256 * (_tmp241 * (-_tmp242 * _tmp267 - _tmp255 * _tmp269) - _tmp243 * _tmp268) +
           _tmp257 * _tmp269 - _tmp268 * p_c(2, 0));
  _res(0, 2) = Scalar(-1.0);
  _res(1, 2) = 0;
  _res(2, 2) = 0;
  _res(3, 2) = 0;

  return _res;
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym

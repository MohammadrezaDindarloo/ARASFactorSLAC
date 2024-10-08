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
 * Symbolic function: IK_residual_func_cost1_wrt_pc_Nl13
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
Eigen::Matrix<Scalar, 4, 3> IkResidualFuncCost1WrtPcNl13(
    const Scalar fh1, const Scalar fv1, const sym::Rot3<Scalar>& DeltaRot,
    const Eigen::Matrix<Scalar, 3, 1>& position_vector, const Eigen::Matrix<Scalar, 4, 1>& encoder,
    const Eigen::Matrix<Scalar, 3, 1>& p_a, const Eigen::Matrix<Scalar, 3, 1>& p_b,
    const Eigen::Matrix<Scalar, 3, 1>& p_c, const Eigen::Matrix<Scalar, 3, 1>& p_d,
    const sym::Rot3<Scalar>& Rot_init, const Scalar epsilon) {
  // Total ops: 776

  // Unused inputs
  (void)encoder;
  (void)epsilon;

  // Input arrays
  const Eigen::Matrix<Scalar, 4, 1>& _DeltaRot = DeltaRot.Data();
  const Eigen::Matrix<Scalar, 4, 1>& _Rot_init = Rot_init.Data();

  // Intermediate terms (272)
  const Scalar _tmp0 = -_DeltaRot[0] * _Rot_init[1] + _DeltaRot[1] * _Rot_init[0] +
                       _DeltaRot[2] * _Rot_init[3] + _DeltaRot[3] * _Rot_init[2];
  const Scalar _tmp1 = -2 * std::pow(_tmp0, Scalar(2));
  const Scalar _tmp2 = _DeltaRot[0] * _Rot_init[2] + _DeltaRot[1] * _Rot_init[3] -
                       _DeltaRot[2] * _Rot_init[0] + _DeltaRot[3] * _Rot_init[1];
  const Scalar _tmp3 = 1 - 2 * std::pow(_tmp2, Scalar(2));
  const Scalar _tmp4 = Scalar(0.20999999999999999) * _tmp1 + Scalar(0.20999999999999999) * _tmp3;
  const Scalar _tmp5 = _DeltaRot[0] * _Rot_init[3] - _DeltaRot[1] * _Rot_init[2] +
                       _DeltaRot[2] * _Rot_init[1] + _DeltaRot[3] * _Rot_init[0];
  const Scalar _tmp6 = 2 * _tmp0 * _tmp5;
  const Scalar _tmp7 = -2 * _DeltaRot[0] * _Rot_init[0] - 2 * _DeltaRot[1] * _Rot_init[1] -
                       2 * _DeltaRot[2] * _Rot_init[2] + 2 * _DeltaRot[3] * _Rot_init[3];
  const Scalar _tmp8 = _tmp2 * _tmp7;
  const Scalar _tmp9 = _tmp6 + _tmp8;
  const Scalar _tmp10 = -Scalar(0.010999999999999999) * _tmp9;
  const Scalar _tmp11 = 2 * _tmp2;
  const Scalar _tmp12 = _tmp11 * _tmp5;
  const Scalar _tmp13 = _tmp0 * _tmp7;
  const Scalar _tmp14 = Scalar(0.20999999999999999) * _tmp12 - Scalar(0.20999999999999999) * _tmp13;
  const Scalar _tmp15 = _tmp10 + _tmp14;
  const Scalar _tmp16 = _tmp15 + _tmp4;
  const Scalar _tmp17 = _tmp16 + position_vector(0, 0);
  const Scalar _tmp18 = -_tmp17 + p_c(0, 0);
  const Scalar _tmp19 = Scalar(1.0) / (fh1);
  const Scalar _tmp20 = Scalar(0.20999999999999999) * _tmp12 + Scalar(0.20999999999999999) * _tmp13;
  const Scalar _tmp21 = -2 * std::pow(_tmp5, Scalar(2));
  const Scalar _tmp22 = Scalar(0.20999999999999999) * _tmp1 + Scalar(0.20999999999999999) * _tmp21 +
                        Scalar(0.20999999999999999);
  const Scalar _tmp23 = _tmp0 * _tmp11;
  const Scalar _tmp24 = _tmp5 * _tmp7;
  const Scalar _tmp25 = _tmp23 - _tmp24;
  const Scalar _tmp26 = -Scalar(0.010999999999999999) * _tmp25;
  const Scalar _tmp27 = _tmp22 + _tmp26;
  const Scalar _tmp28 = _tmp20 + _tmp27;
  const Scalar _tmp29 = _tmp28 + position_vector(1, 0);
  const Scalar _tmp30 = -_tmp29 + p_c(1, 0);
  const Scalar _tmp31 =
      std::sqrt(Scalar(std::pow(_tmp18, Scalar(2)) + std::pow(_tmp30, Scalar(2))));
  const Scalar _tmp32 =
      Scalar(1.0) *
      std::sinh(Scalar(0.71007031138673404) * _tmp19 *
                (-_tmp31 - Scalar(1.4083112389913199) * fh1 * std::asinh(_tmp19 * fv1))) /
      _tmp31;
  const Scalar _tmp33 = _tmp29 - p_c(1, 0);
  const Scalar _tmp34 = Scalar(0.20999999999999999) * _tmp23 + Scalar(0.20999999999999999) * _tmp24;
  const Scalar _tmp35 = -_tmp34;
  const Scalar _tmp36 =
      -Scalar(0.010999999999999999) * _tmp21 - Scalar(0.010999999999999999) * _tmp3;
  const Scalar _tmp37 = Scalar(0.20999999999999999) * _tmp6 - Scalar(0.20999999999999999) * _tmp8;
  const Scalar _tmp38 = _tmp36 + _tmp37;
  const Scalar _tmp39 = _tmp35 + _tmp38;
  const Scalar _tmp40 = _tmp10 - _tmp14;
  const Scalar _tmp41 = _tmp4 + _tmp40;
  const Scalar _tmp42 = _tmp41 + position_vector(0, 0);
  const Scalar _tmp43 = _tmp42 - p_b(0, 0);
  const Scalar _tmp44 = -_tmp22 + _tmp26;
  const Scalar _tmp45 = _tmp20 + _tmp44;
  const Scalar _tmp46 = _tmp45 + position_vector(1, 0);
  const Scalar _tmp47 = _tmp46 - p_b(1, 0);
  const Scalar _tmp48 = std::pow(Scalar(std::pow(_tmp43, Scalar(2)) + std::pow(_tmp47, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp49 = _tmp47 * _tmp48;
  const Scalar _tmp50 = _tmp36 - _tmp37;
  const Scalar _tmp51 = _tmp34 + _tmp50;
  const Scalar _tmp52 = -_tmp4;
  const Scalar _tmp53 = _tmp15 + _tmp52;
  const Scalar _tmp54 = _tmp53 + position_vector(0, 0);
  const Scalar _tmp55 = _tmp54 - p_d(0, 0);
  const Scalar _tmp56 = -_tmp20;
  const Scalar _tmp57 = _tmp27 + _tmp56;
  const Scalar _tmp58 = _tmp57 + position_vector(1, 0);
  const Scalar _tmp59 = _tmp58 - p_d(1, 0);
  const Scalar _tmp60 = std::pow(Scalar(std::pow(_tmp55, Scalar(2)) + std::pow(_tmp59, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp61 = _tmp59 * _tmp60;
  const Scalar _tmp62 = _tmp40 + _tmp52;
  const Scalar _tmp63 = _tmp62 + position_vector(0, 0);
  const Scalar _tmp64 = _tmp63 - p_a(0, 0);
  const Scalar _tmp65 = Scalar(1.0) / (_tmp64);
  const Scalar _tmp66 = _tmp44 + _tmp56;
  const Scalar _tmp67 = _tmp66 + position_vector(1, 0);
  const Scalar _tmp68 = _tmp67 - p_a(1, 0);
  const Scalar _tmp69 = _tmp65 * _tmp68;
  const Scalar _tmp70 = _tmp35 + _tmp50;
  const Scalar _tmp71 = _tmp55 * _tmp60;
  const Scalar _tmp72 = _tmp70 * _tmp71;
  const Scalar _tmp73 = _tmp51 * _tmp61 - _tmp69 * _tmp72;
  const Scalar _tmp74 = Scalar(1.0) / (-_tmp61 + _tmp69 * _tmp71);
  const Scalar _tmp75 = _tmp43 * _tmp48;
  const Scalar _tmp76 = -_tmp49 + _tmp69 * _tmp75;
  const Scalar _tmp77 = _tmp74 * _tmp76;
  const Scalar _tmp78 = _tmp70 * _tmp75;
  const Scalar _tmp79 = _tmp39 * _tmp49 - _tmp69 * _tmp78 - _tmp73 * _tmp77;
  const Scalar _tmp80 = Scalar(1.0) * _tmp66;
  const Scalar _tmp81 = -_tmp80;
  const Scalar _tmp82 = Scalar(1.0) / (_tmp57 + _tmp81);
  const Scalar _tmp83 = Scalar(1.0) * _tmp62;
  const Scalar _tmp84 = -_tmp53 + _tmp83;
  const Scalar _tmp85 = _tmp82 * _tmp84;
  const Scalar _tmp86 = -_tmp51 * _tmp71 + _tmp72;
  const Scalar _tmp87 = -_tmp39 * _tmp75 - _tmp77 * _tmp86 + _tmp78 - _tmp79 * _tmp85;
  const Scalar _tmp88 = Scalar(1.0) / (_tmp87);
  const Scalar _tmp89 = Scalar(1.0) * _tmp82;
  const Scalar _tmp90 = Scalar(1.0) * _tmp74;
  const Scalar _tmp91 = _tmp73 * _tmp74 * _tmp84 * _tmp89 - _tmp86 * _tmp90;
  const Scalar _tmp92 =
      std::sqrt(Scalar(std::pow(_tmp64, Scalar(2)) + std::pow(_tmp68, Scalar(2))));
  const Scalar _tmp93 = Scalar(1.0) / (_tmp92);
  const Scalar _tmp94 = _tmp65 * _tmp92;
  const Scalar _tmp95 = _tmp94 * (_tmp62 * _tmp68 * _tmp93 - _tmp64 * _tmp66 * _tmp93);
  const Scalar _tmp96 = -_tmp53 * _tmp61 + _tmp57 * _tmp71 + _tmp71 * _tmp95;
  const Scalar _tmp97 = -_tmp41 * _tmp49 + _tmp45 * _tmp75 + _tmp75 * _tmp95 - _tmp77 * _tmp96;
  const Scalar _tmp98 = _tmp88 * _tmp97;
  const Scalar _tmp99 = Scalar(1.0) / (_tmp97);
  const Scalar _tmp100 = _tmp87 * _tmp99;
  const Scalar _tmp101 = _tmp100 * (-_tmp90 * _tmp96 - _tmp91 * _tmp98);
  const Scalar _tmp102 = _tmp88 * (_tmp101 + _tmp91);
  const Scalar _tmp103 = _tmp74 * (-_tmp102 * _tmp76 + Scalar(1.0));
  const Scalar _tmp104 = _tmp102 * _tmp75 + _tmp103 * _tmp71;
  const Scalar _tmp105 = _tmp17 - p_c(0, 0);
  const Scalar _tmp106 = std::pow(_tmp105, Scalar(2));
  const Scalar _tmp107 = std::pow(_tmp33, Scalar(2));
  const Scalar _tmp108 = _tmp106 + _tmp107;
  const Scalar _tmp109 = std::pow(_tmp108, Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp110 = _tmp109 * fh1;
  const Scalar _tmp111 = _tmp110 * _tmp94;
  const Scalar _tmp112 = _tmp104 * _tmp111;
  const Scalar _tmp113 = _tmp69 * _tmp74;
  const Scalar _tmp114 = _tmp113 * _tmp73 + _tmp69 * _tmp70;
  const Scalar _tmp115 = _tmp113 * _tmp86 - _tmp114 * _tmp85 - _tmp70;
  const Scalar _tmp116 = _tmp100 * (_tmp113 * _tmp96 - _tmp115 * _tmp98 - _tmp95);
  const Scalar _tmp117 = _tmp115 + _tmp116;
  const Scalar _tmp118 = _tmp117 * _tmp88;
  const Scalar _tmp119 = _tmp74 * (-_tmp118 * _tmp76 - _tmp69);
  const Scalar _tmp120 = _tmp118 * _tmp75 + _tmp119 * _tmp71 + Scalar(1.0);
  const Scalar _tmp121 = _tmp111 * _tmp120;
  const Scalar _tmp122 = Scalar(333.54000000000002) - fv1;
  const Scalar _tmp123 = _tmp80 * _tmp85 + _tmp83;
  const Scalar _tmp124 = 0;
  const Scalar _tmp125 = _tmp124 * _tmp88;
  const Scalar _tmp126 = _tmp109 * _tmp28;
  const Scalar _tmp127 = _tmp109 * _tmp16;
  const Scalar _tmp128 = fh1 * (-_tmp105 * _tmp126 + _tmp127 * _tmp33);
  const Scalar _tmp129 = _tmp76 * _tmp90 * _tmp99;
  const Scalar _tmp130 = Scalar(1.0) * _tmp99;
  const Scalar _tmp131 = _tmp94 * (-_tmp129 * _tmp71 + _tmp130 * _tmp75);
  const Scalar _tmp132 = -_tmp105 * _tmp121 - _tmp112 * _tmp33 -
                         _tmp122 * _tmp94 * (-_tmp125 * _tmp71 * _tmp77 + _tmp125 * _tmp75) -
                         _tmp128 * _tmp131;
  const Scalar _tmp133 = Scalar(1.0) / (_tmp132);
  const Scalar _tmp134 = _tmp45 + _tmp81;
  const Scalar _tmp135 = _tmp134 * _tmp85;
  const Scalar _tmp136 = Scalar(1.0) / (-_tmp135 - _tmp41 + _tmp83);
  const Scalar _tmp137 = Scalar(1.0) * _tmp136;
  const Scalar _tmp138 = _tmp100 * _tmp137;
  const Scalar _tmp139 = -_tmp130 * _tmp79 + _tmp134 * _tmp138;
  const Scalar _tmp140 = Scalar(1.0) * _tmp138 - Scalar(1.0) * _tmp139 * _tmp89;
  const Scalar _tmp141 = _tmp134 * _tmp136;
  const Scalar _tmp142 = _tmp101 * _tmp141 - _tmp102 * _tmp79 - _tmp73 * _tmp90;
  const Scalar _tmp143 = _tmp101 * _tmp137 - _tmp142 * _tmp89;
  const Scalar _tmp144 = Scalar(1.0) * _tmp110;
  const Scalar _tmp145 = _tmp143 * _tmp144;
  const Scalar _tmp146 = fh1 * (_tmp34 + _tmp38);
  const Scalar _tmp147 = _tmp109 * _tmp146;
  const Scalar _tmp148 = -_tmp147 * _tmp33 - Scalar(40.024799999999999) * _tmp25 - _tmp28 * fv1;
  const Scalar _tmp149 = _tmp137 * _tmp85;
  const Scalar _tmp150 = _tmp135 * _tmp137 + Scalar(1.0);
  const Scalar _tmp151 = Scalar(1.0) * _tmp149 - Scalar(1.0) * _tmp150 * _tmp89;
  const Scalar _tmp152 = _tmp105 * _tmp147 + _tmp16 * fv1 + Scalar(40.024799999999999) * _tmp9;
  const Scalar _tmp153 = _tmp134 * _tmp82;
  const Scalar _tmp154 = _tmp137 * _tmp153;
  const Scalar _tmp155 = -Scalar(1.0) * _tmp137 + Scalar(1.0) * _tmp154;
  const Scalar _tmp156 = _tmp79 * _tmp88;
  const Scalar _tmp157 = _tmp114 + _tmp116 * _tmp141 - _tmp117 * _tmp156;
  const Scalar _tmp158 = _tmp116 * _tmp137 - _tmp157 * _tmp89;
  const Scalar _tmp159 = _tmp144 * _tmp158;
  const Scalar _tmp160 = _tmp123 * _tmp136;
  const Scalar _tmp161 = -_tmp124 * _tmp156 - _tmp134 * _tmp160 + _tmp81;
  const Scalar _tmp162 =
      _tmp105 * _tmp159 +
      Scalar(1.0) * _tmp122 * (-_tmp123 * _tmp137 - _tmp161 * _tmp89 + Scalar(1.0)) +
      _tmp128 * _tmp140 + _tmp145 * _tmp33 + _tmp148 * _tmp151 + _tmp152 * _tmp155;
  const Scalar _tmp163 = std::asinh(_tmp133 * _tmp162);
  const Scalar _tmp164 = Scalar(1.0) * _tmp163;
  const Scalar _tmp165 = Scalar(1.0) * std::sinh(_tmp164);
  const Scalar _tmp166 = std::pow(_tmp132, Scalar(-2));
  const Scalar _tmp167 =
      std::pow(Scalar(std::pow(_tmp162, Scalar(2)) * _tmp166 + 1), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp168 = std::pow(_tmp108, Scalar(Scalar(-3) / Scalar(2)));
  const Scalar _tmp169 = _tmp106 * _tmp168;
  const Scalar _tmp170 = _tmp16 * _tmp168;
  const Scalar _tmp171 = _tmp105 * _tmp33;
  const Scalar _tmp172 = _tmp126 - _tmp169 * _tmp28 + _tmp170 * _tmp171;
  const Scalar _tmp173 = _tmp172 * fh1;
  const Scalar _tmp174 = _tmp168 * _tmp171;
  const Scalar _tmp175 = _tmp146 * _tmp174;
  const Scalar _tmp176 = _tmp146 * _tmp169 - _tmp147;
  const Scalar _tmp177 = Scalar(1.0) * _tmp143;
  const Scalar _tmp178 = _tmp174 * fh1;
  const Scalar _tmp179 = _tmp169 * fh1;
  const Scalar _tmp180 = Scalar(1.0) * _tmp158;
  const Scalar _tmp181 = _tmp104 * _tmp94;
  const Scalar _tmp182 = _tmp120 * _tmp94;
  const Scalar _tmp183 = _tmp121 - _tmp131 * _tmp173 - _tmp178 * _tmp181 - _tmp179 * _tmp182;
  const Scalar _tmp184 = _tmp162 * _tmp166;
  const Scalar _tmp185 =
      _tmp167 * (_tmp133 * (_tmp140 * _tmp173 - _tmp151 * _tmp175 + _tmp155 * _tmp176 - _tmp159 +
                            _tmp177 * _tmp178 + _tmp179 * _tmp180) -
                 _tmp183 * _tmp184);
  const Scalar _tmp186 = Scalar(0.71007031138673404) * _tmp166;
  const Scalar _tmp187 = _tmp186 * p_a(2, 0);
  const Scalar _tmp188 = Scalar(1.4083112389913199) * _tmp132;
  const Scalar _tmp189 =
      -_tmp163 * _tmp188 - std::sqrt(Scalar(std::pow(Scalar(-_tmp63 + p_a(0, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp67 + p_a(1, 0)), Scalar(2))));
  const Scalar _tmp190 = Scalar(0.71007031138673404) * _tmp133;
  const Scalar _tmp191 = _tmp189 * _tmp190;
  const Scalar _tmp192 = std::sinh(_tmp191);
  const Scalar _tmp193 = _tmp186 * _tmp189;
  const Scalar _tmp194 = Scalar(1.4083112389913199) * _tmp183;
  const Scalar _tmp195 = _tmp190 * p_a(2, 0) + std::cosh(_tmp164) - std::cosh(_tmp191);
  const Scalar _tmp196 = _tmp110 * _tmp119;
  const Scalar _tmp197 = _tmp103 * _tmp178 + _tmp119 * _tmp179 - _tmp129 * _tmp173 - _tmp196;
  const Scalar _tmp198 = _tmp122 * _tmp125;
  const Scalar _tmp199 = _tmp103 * _tmp110;
  const Scalar _tmp200 =
      _tmp105 * _tmp196 - _tmp128 * _tmp129 - _tmp198 * _tmp77 + _tmp199 * _tmp33;
  const Scalar _tmp201 = Scalar(1.0) / (_tmp200);
  const Scalar _tmp202 = Scalar(0.71007031138673404) * _tmp201;
  const Scalar _tmp203 = _tmp110 * _tmp82;
  const Scalar _tmp204 = _tmp142 * _tmp203;
  const Scalar _tmp205 = _tmp157 * _tmp203;
  const Scalar _tmp206 = _tmp137 * _tmp152;
  const Scalar _tmp207 = _tmp150 * _tmp82;
  const Scalar _tmp208 = _tmp139 * _tmp82;
  const Scalar _tmp209 = _tmp105 * _tmp205 + _tmp122 * _tmp161 * _tmp82 + _tmp128 * _tmp208 +
                         _tmp148 * _tmp207 - _tmp153 * _tmp206 + _tmp204 * _tmp33;
  const Scalar _tmp210 = std::asinh(_tmp201 * _tmp209);
  const Scalar _tmp211 = Scalar(1.4083112389913199) * _tmp210;
  const Scalar _tmp212 =
      -_tmp200 * _tmp211 - std::sqrt(Scalar(std::pow(Scalar(-_tmp54 + p_d(0, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp58 + p_d(1, 0)), Scalar(2))));
  const Scalar _tmp213 = _tmp202 * _tmp212;
  const Scalar _tmp214 = Scalar(1.0) * _tmp210;
  const Scalar _tmp215 = Scalar(1.4083112389913199) * _tmp202 * p_d(2, 0) -
                         Scalar(1.4083112389913199) * std::cosh(_tmp213) +
                         Scalar(1.4083112389913199) * std::cosh(_tmp214);
  const Scalar _tmp216 = std::pow(_tmp200, Scalar(-2));
  const Scalar _tmp217 = _tmp209 * _tmp216;
  const Scalar _tmp218 = _tmp157 * _tmp82;
  const Scalar _tmp219 = _tmp142 * _tmp82;
  const Scalar _tmp220 = _tmp137 * _tmp176;
  const Scalar _tmp221 =
      -_tmp197 * _tmp217 + _tmp201 * (-_tmp153 * _tmp220 + _tmp173 * _tmp208 - _tmp175 * _tmp207 +
                                      _tmp178 * _tmp219 + _tmp179 * _tmp218 - _tmp205);
  const Scalar _tmp222 =
      std::pow(Scalar(std::pow(_tmp209, Scalar(2)) * _tmp216 + 1), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp223 = Scalar(1.0) * _tmp222 * std::sinh(_tmp214);
  const Scalar _tmp224 = Scalar(0.71007031138673404) * _tmp216;
  const Scalar _tmp225 = _tmp197 * _tmp224;
  const Scalar _tmp226 = std::sinh(_tmp213);
  const Scalar _tmp227 = Scalar(1.4083112389913199) * _tmp200;
  const Scalar _tmp228 = _tmp222 * _tmp227;
  const Scalar _tmp229 = _tmp110 * _tmp136;
  const Scalar _tmp230 = _tmp101 * _tmp229;
  const Scalar _tmp231 = _tmp116 * _tmp229;
  const Scalar _tmp232 = -_tmp105 * _tmp231 + _tmp122 * _tmp160 - _tmp128 * _tmp138 -
                         _tmp148 * _tmp149 + _tmp206 - _tmp230 * _tmp33;
  const Scalar _tmp233 = _tmp102 * _tmp110;
  const Scalar _tmp234 = _tmp110 * _tmp118;
  const Scalar _tmp235 = _tmp105 * _tmp234 + _tmp128 * _tmp130 + _tmp198 + _tmp233 * _tmp33;
  const Scalar _tmp236 = Scalar(1.0) / (_tmp235);
  const Scalar _tmp237 = std::asinh(_tmp232 * _tmp236);
  const Scalar _tmp238 = Scalar(1.0) * _tmp237;
  const Scalar _tmp239 = Scalar(1.4083112389913199) * _tmp235;
  const Scalar _tmp240 =
      -_tmp237 * _tmp239 - std::sqrt(Scalar(std::pow(Scalar(-_tmp42 + p_b(0, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp46 + p_b(1, 0)), Scalar(2))));
  const Scalar _tmp241 = Scalar(0.71007031138673404) * _tmp236;
  const Scalar _tmp242 = _tmp240 * _tmp241;
  const Scalar _tmp243 = _tmp241 * p_b(2, 0) + std::cosh(_tmp238) - std::cosh(_tmp242);
  const Scalar _tmp244 = _tmp130 * fh1;
  const Scalar _tmp245 = _tmp102 * _tmp178 + _tmp118 * _tmp179 + _tmp172 * _tmp244 - _tmp234;
  const Scalar _tmp246 = Scalar(1.4083112389913199) * _tmp245;
  const Scalar _tmp247 = std::sinh(_tmp242);
  const Scalar _tmp248 = std::pow(_tmp235, Scalar(-2));
  const Scalar _tmp249 =
      std::pow(Scalar(std::pow(_tmp232, Scalar(2)) * _tmp248 + 1), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp250 = _tmp232 * _tmp248;
  const Scalar _tmp251 = _tmp116 * _tmp136;
  const Scalar _tmp252 = _tmp101 * _tmp136;
  const Scalar _tmp253 =
      _tmp249 * (_tmp236 * (-_tmp138 * _tmp173 + _tmp149 * _tmp175 - _tmp178 * _tmp252 -
                            _tmp179 * _tmp251 + _tmp220 + _tmp231) -
                 _tmp245 * _tmp250);
  const Scalar _tmp254 = Scalar(0.71007031138673404) * _tmp248;
  const Scalar _tmp255 = _tmp245 * _tmp254;
  const Scalar _tmp256 = Scalar(1.0) * std::sinh(_tmp238);
  const Scalar _tmp257 = _tmp107 * _tmp168;
  const Scalar _tmp258 = _tmp257 * fh1;
  const Scalar _tmp259 = _tmp107 * _tmp170 - _tmp127 - _tmp174 * _tmp28;
  const Scalar _tmp260 = _tmp259 * fh1;
  const Scalar _tmp261 = _tmp112 - _tmp131 * _tmp260 - _tmp178 * _tmp182 - _tmp181 * _tmp258;
  const Scalar _tmp262 = Scalar(1.4083112389913199) * _tmp261;
  const Scalar _tmp263 = -_tmp146 * _tmp257 + _tmp147;
  const Scalar _tmp264 =
      _tmp167 * (_tmp133 * (_tmp140 * _tmp260 - _tmp145 + _tmp151 * _tmp263 + _tmp155 * _tmp175 +
                            _tmp177 * _tmp258 + _tmp178 * _tmp180) -
                 _tmp184 * _tmp261);
  const Scalar _tmp265 = _tmp103 * _tmp258 + _tmp119 * _tmp178 - _tmp129 * _tmp260 - _tmp199;
  const Scalar _tmp266 = _tmp201 * (-_tmp154 * _tmp175 + _tmp178 * _tmp218 - _tmp204 +
                                    _tmp207 * _tmp263 + _tmp208 * _tmp260 + _tmp219 * _tmp258) -
                         _tmp217 * _tmp265;
  const Scalar _tmp267 = _tmp224 * _tmp265;
  const Scalar _tmp268 = _tmp102 * _tmp258 + _tmp118 * _tmp178 - _tmp233 + _tmp244 * _tmp259;
  const Scalar _tmp269 =
      _tmp249 * (_tmp236 * (_tmp137 * _tmp175 - _tmp138 * _tmp260 - _tmp149 * _tmp263 -
                            _tmp178 * _tmp251 + _tmp230 - _tmp252 * _tmp258) -
                 _tmp250 * _tmp268);
  const Scalar _tmp270 = _tmp254 * _tmp268;
  const Scalar _tmp271 = Scalar(1.4083112389913199) * _tmp268;

  // Output terms (1)
  Eigen::Matrix<Scalar, 4, 3> _res;

  _res(0, 0) = -_tmp18 * _tmp32;
  _res(1, 0) =
      -_tmp188 *
          (_tmp165 * _tmp185 - _tmp183 * _tmp187 -
           _tmp192 * (-_tmp183 * _tmp193 + _tmp190 * (-_tmp163 * _tmp194 - _tmp185 * _tmp188))) -
      _tmp194 * _tmp195;
  _res(2, 0) =
      -_tmp197 * _tmp215 -
      _tmp227 *
          (_tmp221 * _tmp223 - _tmp225 * p_d(2, 0) -
           _tmp226 * (_tmp202 * (-_tmp197 * _tmp211 - _tmp221 * _tmp228) - _tmp212 * _tmp225));
  _res(3, 0) =
      -_tmp239 *
          (-_tmp247 * (-_tmp240 * _tmp255 + _tmp241 * (-_tmp237 * _tmp246 - _tmp239 * _tmp253)) +
           _tmp253 * _tmp256 - _tmp255 * p_b(2, 0)) -
      _tmp243 * _tmp246;
  _res(0, 1) = -_tmp30 * _tmp32;
  _res(1, 1) =
      -_tmp188 *
          (_tmp165 * _tmp264 - _tmp187 * _tmp261 -
           _tmp192 * (_tmp190 * (-_tmp163 * _tmp262 - _tmp188 * _tmp264) - _tmp193 * _tmp261)) -
      _tmp195 * _tmp262;
  _res(2, 1) =
      -_tmp215 * _tmp265 -
      _tmp227 *
          (_tmp223 * _tmp266 -
           _tmp226 * (_tmp202 * (-_tmp211 * _tmp265 - _tmp228 * _tmp266) - _tmp212 * _tmp267) -
           _tmp267 * p_d(2, 0));
  _res(3, 1) =
      -_tmp239 *
          (-_tmp247 * (-_tmp240 * _tmp270 + _tmp241 * (-_tmp237 * _tmp271 - _tmp239 * _tmp269)) +
           _tmp256 * _tmp269 - _tmp270 * p_b(2, 0)) -
      _tmp243 * _tmp271;
  _res(0, 2) = Scalar(-1.0);
  _res(1, 2) = 0;
  _res(2, 2) = 0;
  _res(3, 2) = 0;

  return _res;
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym

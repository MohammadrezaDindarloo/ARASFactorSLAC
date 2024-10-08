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
 * Symbolic function: IK_residual_func_cost1_wrt_pa_Nl11
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
Eigen::Matrix<Scalar, 4, 3> IkResidualFuncCost1WrtPaNl11(
    const Scalar fh1, const Scalar fv1, const sym::Rot3<Scalar>& DeltaRot,
    const Eigen::Matrix<Scalar, 3, 1>& position_vector, const Eigen::Matrix<Scalar, 4, 1>& encoder,
    const Eigen::Matrix<Scalar, 3, 1>& p_a, const Eigen::Matrix<Scalar, 3, 1>& p_b,
    const Eigen::Matrix<Scalar, 3, 1>& p_c, const Eigen::Matrix<Scalar, 3, 1>& p_d,
    const sym::Rot3<Scalar>& Rot_init, const Scalar epsilon) {
  // Total ops: 1058

  // Unused inputs
  (void)encoder;
  (void)epsilon;

  // Input arrays
  const Eigen::Matrix<Scalar, 4, 1>& _DeltaRot = DeltaRot.Data();
  const Eigen::Matrix<Scalar, 4, 1>& _Rot_init = Rot_init.Data();

  // Intermediate terms (329)
  const Scalar _tmp0 = _DeltaRot[0] * _Rot_init[2] + _DeltaRot[1] * _Rot_init[3] -
                       _DeltaRot[2] * _Rot_init[0] + _DeltaRot[3] * _Rot_init[1];
  const Scalar _tmp1 = -2 * std::pow(_tmp0, Scalar(2));
  const Scalar _tmp2 = -_DeltaRot[0] * _Rot_init[1] + _DeltaRot[1] * _Rot_init[0] +
                       _DeltaRot[2] * _Rot_init[3] + _DeltaRot[3] * _Rot_init[2];
  const Scalar _tmp3 = 1 - 2 * std::pow(_tmp2, Scalar(2));
  const Scalar _tmp4 = Scalar(0.20999999999999999) * _tmp1 + Scalar(0.20999999999999999) * _tmp3;
  const Scalar _tmp5 = -_tmp4;
  const Scalar _tmp6 = _DeltaRot[0] * _Rot_init[3] - _DeltaRot[1] * _Rot_init[2] +
                       _DeltaRot[2] * _Rot_init[1] + _DeltaRot[3] * _Rot_init[0];
  const Scalar _tmp7 = 2 * _tmp6;
  const Scalar _tmp8 = _tmp2 * _tmp7;
  const Scalar _tmp9 = -2 * _DeltaRot[0] * _Rot_init[0] - 2 * _DeltaRot[1] * _Rot_init[1] -
                       2 * _DeltaRot[2] * _Rot_init[2] + 2 * _DeltaRot[3] * _Rot_init[3];
  const Scalar _tmp10 = _tmp0 * _tmp9;
  const Scalar _tmp11 = _tmp10 + _tmp8;
  const Scalar _tmp12 = -Scalar(0.010999999999999999) * _tmp11;
  const Scalar _tmp13 = _tmp0 * _tmp7;
  const Scalar _tmp14 = _tmp2 * _tmp9;
  const Scalar _tmp15 = Scalar(0.20999999999999999) * _tmp13 - Scalar(0.20999999999999999) * _tmp14;
  const Scalar _tmp16 = _tmp12 - _tmp15;
  const Scalar _tmp17 = _tmp16 + _tmp5;
  const Scalar _tmp18 = _tmp17 + position_vector(0, 0);
  const Scalar _tmp19 = _tmp18 - p_a(0, 0);
  const Scalar _tmp20 = std::pow(_tmp19, Scalar(2));
  const Scalar _tmp21 = -2 * std::pow(_tmp6, Scalar(2));
  const Scalar _tmp22 = Scalar(0.20999999999999999) * _tmp21 + Scalar(0.20999999999999999) * _tmp3;
  const Scalar _tmp23 = -_tmp22;
  const Scalar _tmp24 = Scalar(0.20999999999999999) * _tmp13 + Scalar(0.20999999999999999) * _tmp14;
  const Scalar _tmp25 = 2 * _tmp0 * _tmp2;
  const Scalar _tmp26 = _tmp6 * _tmp9;
  const Scalar _tmp27 = _tmp25 - _tmp26;
  const Scalar _tmp28 = -Scalar(0.010999999999999999) * _tmp27;
  const Scalar _tmp29 = -_tmp24 + _tmp28;
  const Scalar _tmp30 = _tmp23 + _tmp29;
  const Scalar _tmp31 = _tmp30 + position_vector(1, 0);
  const Scalar _tmp32 = _tmp31 - p_a(1, 0);
  const Scalar _tmp33 = std::pow(_tmp32, Scalar(2));
  const Scalar _tmp34 = _tmp20 + _tmp33;
  const Scalar _tmp35 = std::pow(_tmp34, Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp36 = -Scalar(0.20999999999999999) * _tmp10 + Scalar(0.20999999999999999) * _tmp8;
  const Scalar _tmp37 = -_tmp36;
  const Scalar _tmp38 = -Scalar(0.010999999999999999) * _tmp1 -
                        Scalar(0.010999999999999999) * _tmp21 + Scalar(-0.010999999999999999);
  const Scalar _tmp39 = Scalar(0.20999999999999999) * _tmp25 + Scalar(0.20999999999999999) * _tmp26;
  const Scalar _tmp40 = _tmp38 - _tmp39;
  const Scalar _tmp41 = _tmp37 + _tmp40;
  const Scalar _tmp42 = _tmp35 * _tmp41;
  const Scalar _tmp43 = _tmp12 + _tmp15;
  const Scalar _tmp44 = _tmp43 + _tmp5;
  const Scalar _tmp45 = _tmp44 + position_vector(0, 0);
  const Scalar _tmp46 = _tmp45 - p_d(0, 0);
  const Scalar _tmp47 = Scalar(1.0) / (_tmp46);
  const Scalar _tmp48 = _tmp22 + _tmp29;
  const Scalar _tmp49 = _tmp48 + position_vector(1, 0);
  const Scalar _tmp50 = _tmp49 - p_d(1, 0);
  const Scalar _tmp51 = _tmp47 * _tmp50;
  const Scalar _tmp52 = _tmp35 * _tmp51;
  const Scalar _tmp53 = _tmp32 * _tmp35;
  const Scalar _tmp54 = _tmp19 * _tmp52 - _tmp53;
  const Scalar _tmp55 = _tmp4 + _tmp43;
  const Scalar _tmp56 = _tmp55 + position_vector(0, 0);
  const Scalar _tmp57 = _tmp56 - p_c(0, 0);
  const Scalar _tmp58 = _tmp24 + _tmp28;
  const Scalar _tmp59 = _tmp22 + _tmp58;
  const Scalar _tmp60 = _tmp59 + position_vector(1, 0);
  const Scalar _tmp61 = _tmp60 - p_c(1, 0);
  const Scalar _tmp62 = std::pow(Scalar(std::pow(_tmp57, Scalar(2)) + std::pow(_tmp61, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp63 = _tmp57 * _tmp62;
  const Scalar _tmp64 = _tmp61 * _tmp62;
  const Scalar _tmp65 = Scalar(1.0) / (_tmp51 * _tmp63 - _tmp64);
  const Scalar _tmp66 = _tmp38 + _tmp39;
  const Scalar _tmp67 = _tmp37 + _tmp66;
  const Scalar _tmp68 = _tmp63 * _tmp67;
  const Scalar _tmp69 = _tmp36 + _tmp66;
  const Scalar _tmp70 = _tmp65 * (-_tmp51 * _tmp68 + _tmp64 * _tmp69);
  const Scalar _tmp71 = _tmp35 * _tmp67;
  const Scalar _tmp72 = _tmp19 * _tmp71;
  const Scalar _tmp73 = _tmp41 * _tmp53 - _tmp51 * _tmp72 - _tmp54 * _tmp70;
  const Scalar _tmp74 = Scalar(1.0) * _tmp48;
  const Scalar _tmp75 = -_tmp74;
  const Scalar _tmp76 = Scalar(1.0) / (_tmp59 + _tmp75);
  const Scalar _tmp77 = Scalar(1.0) * _tmp44;
  const Scalar _tmp78 = _tmp76 * (-_tmp55 + _tmp77);
  const Scalar _tmp79 = _tmp65 * (-_tmp63 * _tmp69 + _tmp68);
  const Scalar _tmp80 = -_tmp19 * _tmp42 - _tmp54 * _tmp79 + _tmp72 - _tmp73 * _tmp78;
  const Scalar _tmp81 = Scalar(1.0) / (_tmp80);
  const Scalar _tmp82 = _tmp74 * _tmp78 + _tmp77;
  const Scalar _tmp83 = 0;
  const Scalar _tmp84 = _tmp81 * _tmp83;
  const Scalar _tmp85 = _tmp35 * _tmp84;
  const Scalar _tmp86 = std::pow(_tmp34, Scalar(Scalar(-3) / Scalar(2)));
  const Scalar _tmp87 = _tmp20 * _tmp86;
  const Scalar _tmp88 = _tmp19 * _tmp32 * _tmp86;
  const Scalar _tmp89 = _tmp51 * _tmp87 - _tmp52 - _tmp88;
  const Scalar _tmp90 = _tmp67 * _tmp87;
  const Scalar _tmp91 = _tmp41 * _tmp88;
  const Scalar _tmp92 = _tmp51 * _tmp71 - _tmp51 * _tmp90 - _tmp70 * _tmp89 + _tmp91;
  const Scalar _tmp93 =
      -_tmp41 * _tmp87 + _tmp42 - _tmp71 - _tmp78 * _tmp92 - _tmp79 * _tmp89 + _tmp90;
  const Scalar _tmp94 = std::pow(_tmp80, Scalar(-2));
  const Scalar _tmp95 = _tmp93 * _tmp94;
  const Scalar _tmp96 = _tmp19 * _tmp35;
  const Scalar _tmp97 = _tmp83 * _tmp96;
  const Scalar _tmp98 = _tmp63 * _tmp65;
  const Scalar _tmp99 = _tmp84 * _tmp98;
  const Scalar _tmp100 = _tmp83 * _tmp98;
  const Scalar _tmp101 = _tmp100 * _tmp54;
  const Scalar _tmp102 = Scalar(333.54000000000002) - fv1;
  const Scalar _tmp103 =
      std::sqrt(Scalar(std::pow(_tmp46, Scalar(2)) + std::pow(_tmp50, Scalar(2))));
  const Scalar _tmp104 = _tmp103 * _tmp47;
  const Scalar _tmp105 = _tmp102 * _tmp104;
  const Scalar _tmp106 = Scalar(1.0) / (_tmp103);
  const Scalar _tmp107 = _tmp104 * (_tmp106 * _tmp44 * _tmp50 - _tmp106 * _tmp46 * _tmp48);
  const Scalar _tmp108 = _tmp65 * (_tmp107 * _tmp63 - _tmp55 * _tmp64 + _tmp59 * _tmp63);
  const Scalar _tmp109 = _tmp30 * _tmp35;
  const Scalar _tmp110 = _tmp107 * _tmp96 - _tmp108 * _tmp54 + _tmp109 * _tmp19 - _tmp17 * _tmp53;
  const Scalar _tmp111 = _tmp51 * _tmp67;
  const Scalar _tmp112 = _tmp111 + _tmp51 * _tmp70;
  const Scalar _tmp113 = -_tmp112 * _tmp78 + _tmp51 * _tmp79 - _tmp67;
  const Scalar _tmp114 = _tmp113 * _tmp81;
  const Scalar _tmp115 = -_tmp107 + _tmp108 * _tmp51 - _tmp110 * _tmp114;
  const Scalar _tmp116 = Scalar(1.0) / (_tmp110);
  const Scalar _tmp117 = _tmp116 * _tmp80;
  const Scalar _tmp118 = _tmp115 * _tmp117;
  const Scalar _tmp119 = _tmp113 + _tmp118;
  const Scalar _tmp120 = _tmp119 * _tmp81;
  const Scalar _tmp121 = _tmp120 * _tmp35;
  const Scalar _tmp122 = _tmp115 * _tmp116;
  const Scalar _tmp123 = _tmp122 * _tmp93;
  const Scalar _tmp124 = -_tmp107 * _tmp35 + _tmp107 * _tmp87 - _tmp108 * _tmp89 - _tmp109 -
                         _tmp17 * _tmp88 + _tmp30 * _tmp87;
  const Scalar _tmp125 = _tmp110 * _tmp113;
  const Scalar _tmp126 = _tmp117 * (-_tmp114 * _tmp124 + _tmp125 * _tmp95);
  const Scalar _tmp127 = _tmp115 * _tmp80;
  const Scalar _tmp128 = std::pow(_tmp110, Scalar(-2));
  const Scalar _tmp129 = _tmp124 * _tmp128;
  const Scalar _tmp130 = _tmp127 * _tmp129;
  const Scalar _tmp131 = _tmp123 + _tmp126 - _tmp130;
  const Scalar _tmp132 = _tmp131 * _tmp81;
  const Scalar _tmp133 = _tmp119 * _tmp95;
  const Scalar _tmp134 = _tmp54 * _tmp81;
  const Scalar _tmp135 = -_tmp120 * _tmp89 - _tmp131 * _tmp134 + _tmp133 * _tmp54;
  const Scalar _tmp136 = _tmp16 + _tmp4;
  const Scalar _tmp137 = _tmp136 - p_b(0, 0) + position_vector(0, 0);
  const Scalar _tmp138 = _tmp23 + _tmp58;
  const Scalar _tmp139 = _tmp138 - p_b(1, 0) + position_vector(1, 0);
  const Scalar _tmp140 =
      std::pow(Scalar(std::pow(_tmp137, Scalar(2)) + std::pow(_tmp139, Scalar(2))),
               Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp141 = _tmp137 * _tmp140;
  const Scalar _tmp142 = _tmp141 * fh1;
  const Scalar _tmp143 = _tmp104 * _tmp142;
  const Scalar _tmp144 = Scalar(1.0) * _tmp70;
  const Scalar _tmp145 = _tmp144 * _tmp78 - Scalar(1.0) * _tmp79;
  const Scalar _tmp146 = _tmp145 * _tmp81;
  const Scalar _tmp147 = -Scalar(1.0) * _tmp108 - _tmp110 * _tmp146;
  const Scalar _tmp148 = _tmp116 * _tmp147;
  const Scalar _tmp149 = _tmp148 * _tmp93;
  const Scalar _tmp150 = _tmp110 * _tmp145;
  const Scalar _tmp151 = _tmp117 * (-_tmp124 * _tmp146 + _tmp150 * _tmp95);
  const Scalar _tmp152 = _tmp147 * _tmp80;
  const Scalar _tmp153 = _tmp129 * _tmp152;
  const Scalar _tmp154 = _tmp149 + _tmp151 - _tmp153;
  const Scalar _tmp155 = _tmp81 * _tmp96;
  const Scalar _tmp156 = _tmp117 * _tmp147;
  const Scalar _tmp157 = _tmp145 + _tmp156;
  const Scalar _tmp158 = _tmp157 * _tmp95;
  const Scalar _tmp159 = _tmp157 * _tmp81;
  const Scalar _tmp160 = -_tmp134 * _tmp154 + _tmp158 * _tmp54 - _tmp159 * _tmp89;
  const Scalar _tmp161 = _tmp159 * _tmp35;
  const Scalar _tmp162 = _tmp139 * _tmp140;
  const Scalar _tmp163 = _tmp162 * fh1;
  const Scalar _tmp164 = _tmp104 * _tmp163;
  const Scalar _tmp165 = Scalar(1.0) * _tmp116;
  const Scalar _tmp166 = _tmp165 * _tmp98;
  const Scalar _tmp167 = Scalar(1.0) * _tmp54 * _tmp98;
  const Scalar _tmp168 = Scalar(1.0) * _tmp96;
  const Scalar _tmp169 = _tmp165 * _tmp35;
  const Scalar _tmp170 = fh1 * (_tmp136 * _tmp162 - _tmp138 * _tmp141);
  const Scalar _tmp171 = _tmp104 * _tmp170;
  const Scalar _tmp172 =
      -_tmp105 * (_tmp101 * _tmp95 + _tmp84 * _tmp87 - _tmp85 - _tmp89 * _tmp99 - _tmp95 * _tmp97) -
      _tmp143 *
          (_tmp120 * _tmp87 - _tmp121 + _tmp132 * _tmp96 - _tmp133 * _tmp96 + _tmp135 * _tmp98) -
      _tmp164 *
          (_tmp154 * _tmp155 - _tmp158 * _tmp96 + _tmp159 * _tmp87 + _tmp160 * _tmp98 - _tmp161) -
      _tmp171 *
          (_tmp129 * _tmp167 - _tmp129 * _tmp168 + _tmp165 * _tmp87 - _tmp166 * _tmp89 - _tmp169);
  const Scalar _tmp173 = _tmp30 + _tmp75;
  const Scalar _tmp174 = _tmp173 * _tmp78;
  const Scalar _tmp175 = Scalar(1.0) / (-_tmp17 - _tmp174 + _tmp77);
  const Scalar _tmp176 = Scalar(1.0) * _tmp175;
  const Scalar _tmp177 = _tmp173 * _tmp176;
  const Scalar _tmp178 = fh1 * (_tmp36 + _tmp40);
  const Scalar _tmp179 = Scalar(40.024799999999999) * _tmp11 + _tmp136 * fv1 + _tmp141 * _tmp178;
  const Scalar _tmp180 = _tmp173 * _tmp175;
  const Scalar _tmp181 = _tmp73 * _tmp81;
  const Scalar _tmp182 = _tmp112 + _tmp118 * _tmp180 - _tmp119 * _tmp181;
  const Scalar _tmp183 = Scalar(1.0) * _tmp76;
  const Scalar _tmp184 = Scalar(1.0) * _tmp142;
  const Scalar _tmp185 = _tmp175 * _tmp82;
  const Scalar _tmp186 = -_tmp173 * _tmp185 - _tmp181 * _tmp83 + _tmp75;
  const Scalar _tmp187 = -_tmp144 + _tmp156 * _tmp180 - _tmp157 * _tmp181;
  const Scalar _tmp188 = Scalar(1.0) * _tmp163;
  const Scalar _tmp189 = -_tmp138 * fv1 - _tmp162 * _tmp178 - Scalar(40.024799999999999) * _tmp27;
  const Scalar _tmp190 = _tmp176 * _tmp78;
  const Scalar _tmp191 = _tmp76 * (_tmp174 * _tmp176 + Scalar(1.0));
  const Scalar _tmp192 = _tmp117 * _tmp176;
  const Scalar _tmp193 = -_tmp165 * _tmp73 + _tmp173 * _tmp192;
  const Scalar _tmp194 = Scalar(1.0) * _tmp170;
  const Scalar _tmp195 =
      Scalar(1.0) * _tmp102 * (-_tmp176 * _tmp82 - _tmp183 * _tmp186 + Scalar(1.0)) +
      Scalar(1.0) * _tmp179 * (-_tmp176 + _tmp177 * _tmp76) +
      _tmp184 * (_tmp118 * _tmp176 - _tmp182 * _tmp183) +
      _tmp188 * (_tmp156 * _tmp176 - _tmp183 * _tmp187) +
      Scalar(1.0) * _tmp189 * (_tmp190 - Scalar(1.0) * _tmp191) +
      _tmp194 * (-_tmp183 * _tmp193 + _tmp192);
  const Scalar _tmp196 = -_tmp134 * _tmp157 + Scalar(1.0);
  const Scalar _tmp197 = -_tmp119 * _tmp134 - _tmp51;
  const Scalar _tmp198 = -_tmp105 * (-_tmp100 * _tmp134 + _tmp19 * _tmp85) -
                         _tmp143 * (_tmp121 * _tmp19 + _tmp197 * _tmp98 + Scalar(1.0)) -
                         _tmp164 * (_tmp161 * _tmp19 + _tmp196 * _tmp98) -
                         _tmp171 * (-_tmp166 * _tmp54 + _tmp169 * _tmp19);
  const Scalar _tmp199 = std::pow(_tmp198, Scalar(-2));
  const Scalar _tmp200 = _tmp195 * _tmp199;
  const Scalar _tmp201 = Scalar(1.0) / (_tmp198);
  const Scalar _tmp202 = _tmp149 * _tmp180 + _tmp151 * _tmp180 - _tmp153 * _tmp180 -
                         _tmp154 * _tmp181 + _tmp158 * _tmp73 - _tmp159 * _tmp92;
  const Scalar _tmp203 = _tmp177 * _tmp80;
  const Scalar _tmp204 = Scalar(1.0) * _tmp73;
  const Scalar _tmp205 = _tmp116 * _tmp177;
  const Scalar _tmp206 =
      -_tmp129 * _tmp203 + _tmp129 * _tmp204 - _tmp165 * _tmp92 + _tmp205 * _tmp93;
  const Scalar _tmp207 = _tmp176 * _tmp80;
  const Scalar _tmp208 = _tmp129 * _tmp207;
  const Scalar _tmp209 = _tmp116 * _tmp176;
  const Scalar _tmp210 = _tmp209 * _tmp93;
  const Scalar _tmp211 = _tmp73 * _tmp83;
  const Scalar _tmp212 = _tmp102 * _tmp76;
  const Scalar _tmp213 = _tmp212 * (_tmp211 * _tmp95 - _tmp84 * _tmp92);
  const Scalar _tmp214 = -_tmp120 * _tmp92 + _tmp123 * _tmp180 + _tmp126 * _tmp180 -
                         _tmp130 * _tmp180 - _tmp131 * _tmp181 + _tmp133 * _tmp73;
  const Scalar _tmp215 =
      -_tmp172 * _tmp200 +
      _tmp201 * (_tmp184 * (_tmp123 * _tmp176 + _tmp126 * _tmp176 - _tmp130 * _tmp176 -
                            _tmp183 * _tmp214) +
                 _tmp188 * (_tmp149 * _tmp176 + _tmp151 * _tmp176 - _tmp153 * _tmp176 -
                            _tmp183 * _tmp202) +
                 _tmp194 * (-_tmp183 * _tmp206 - _tmp208 + _tmp210) - Scalar(1.0) * _tmp213);
  const Scalar _tmp216 =
      std::pow(Scalar(std::pow(_tmp195, Scalar(2)) * _tmp199 + 1), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp217 = std::asinh(_tmp195 * _tmp201);
  const Scalar _tmp218 = Scalar(1.0) * _tmp217;
  const Scalar _tmp219 = Scalar(1.0) * _tmp216 * std::sinh(_tmp218);
  const Scalar _tmp220 = Scalar(0.71007031138673404) * _tmp199;
  const Scalar _tmp221 = _tmp220 * p_d(2, 0);
  const Scalar _tmp222 = Scalar(1.4083112389913199) * _tmp198;
  const Scalar _tmp223 =
      -_tmp217 * _tmp222 - std::sqrt(Scalar(std::pow(Scalar(-_tmp45 + p_d(0, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp49 + p_d(1, 0)), Scalar(2))));
  const Scalar _tmp224 = Scalar(0.71007031138673404) * _tmp201;
  const Scalar _tmp225 = _tmp223 * _tmp224;
  const Scalar _tmp226 = std::sinh(_tmp225);
  const Scalar _tmp227 = _tmp220 * _tmp223;
  const Scalar _tmp228 = _tmp216 * _tmp222;
  const Scalar _tmp229 = Scalar(1.4083112389913199) * _tmp172;
  const Scalar _tmp230 = _tmp224 * p_d(2, 0) + std::cosh(_tmp218) - std::cosh(_tmp225);
  const Scalar _tmp231 = _tmp102 * _tmp84;
  const Scalar _tmp232 = _tmp231 * _tmp65;
  const Scalar _tmp233 = _tmp102 * _tmp83;
  const Scalar _tmp234 = _tmp233 * _tmp95;
  const Scalar _tmp235 = _tmp54 * _tmp65;
  const Scalar _tmp236 = _tmp165 * _tmp170;
  const Scalar _tmp237 = _tmp236 * _tmp65;
  const Scalar _tmp238 = _tmp129 * _tmp194;
  const Scalar _tmp239 = _tmp163 * _tmp65;
  const Scalar _tmp240 = _tmp142 * _tmp65;
  const Scalar _tmp241 = _tmp135 * _tmp240 + _tmp160 * _tmp239 - _tmp232 * _tmp89 +
                         _tmp234 * _tmp235 + _tmp235 * _tmp238 - _tmp237 * _tmp89;
  const Scalar _tmp242 =
      -_tmp134 * _tmp233 * _tmp65 + _tmp196 * _tmp239 + _tmp197 * _tmp240 - _tmp235 * _tmp236;
  const Scalar _tmp243 = std::pow(_tmp242, Scalar(-2));
  const Scalar _tmp244 = Scalar(0.71007031138673404) * _tmp243;
  const Scalar _tmp245 = _tmp241 * _tmp244;
  const Scalar _tmp246 = _tmp176 * _tmp179;
  const Scalar _tmp247 = _tmp170 * _tmp76;
  const Scalar _tmp248 = _tmp163 * _tmp76;
  const Scalar _tmp249 = _tmp142 * _tmp76;
  const Scalar _tmp250 = -_tmp173 * _tmp246 * _tmp76 + _tmp182 * _tmp249 + _tmp186 * _tmp212 +
                         _tmp187 * _tmp248 + _tmp189 * _tmp191 + _tmp193 * _tmp247;
  const Scalar _tmp251 = Scalar(1.0) / (_tmp242);
  const Scalar _tmp252 = std::asinh(_tmp250 * _tmp251);
  const Scalar _tmp253 = Scalar(1.0) * _tmp252;
  const Scalar _tmp254 = Scalar(1.0) * std::sinh(_tmp253);
  const Scalar _tmp255 =
      std::pow(Scalar(_tmp243 * std::pow(_tmp250, Scalar(2)) + 1), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp256 = _tmp243 * _tmp250;
  const Scalar _tmp257 =
      _tmp255 * (-_tmp241 * _tmp256 +
                 _tmp251 * (_tmp202 * _tmp248 + _tmp206 * _tmp247 + _tmp213 + _tmp214 * _tmp249));
  const Scalar _tmp258 = Scalar(1.4083112389913199) * _tmp252;
  const Scalar _tmp259 =
      -_tmp242 * _tmp258 - std::sqrt(Scalar(std::pow(Scalar(-_tmp56 + p_c(0, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp60 + p_c(1, 0)), Scalar(2))));
  const Scalar _tmp260 = Scalar(0.71007031138673404) * _tmp251;
  const Scalar _tmp261 = _tmp259 * _tmp260;
  const Scalar _tmp262 = std::sinh(_tmp261);
  const Scalar _tmp263 = Scalar(1.4083112389913199) * _tmp242;
  const Scalar _tmp264 = Scalar(1.4083112389913199) * _tmp260 * p_c(2, 0) +
                         Scalar(1.4083112389913199) * std::cosh(_tmp253) -
                         Scalar(1.4083112389913199) * std::cosh(_tmp261);
  const Scalar _tmp265 = -_tmp31 + p_a(1, 0);
  const Scalar _tmp266 = -_tmp18 + p_a(0, 0);
  const Scalar _tmp267 =
      std::sqrt(Scalar(std::pow(_tmp265, Scalar(2)) + std::pow(_tmp266, Scalar(2))));
  const Scalar _tmp268 = _tmp120 * _tmp142 + _tmp159 * _tmp163 + _tmp231 + _tmp236;
  const Scalar _tmp269 = _tmp163 * _tmp175;
  const Scalar _tmp270 = _tmp142 * _tmp175;
  const Scalar _tmp271 = _tmp102 * _tmp185 - _tmp118 * _tmp270 - _tmp156 * _tmp269 -
                         _tmp170 * _tmp192 - _tmp189 * _tmp190 + _tmp246;
  const Scalar _tmp272 = Scalar(1.0) / (_tmp268);
  const Scalar _tmp273 = std::asinh(_tmp271 * _tmp272);
  const Scalar _tmp274 = Scalar(1.4083112389913199) * _tmp273;
  const Scalar _tmp275 = -_tmp267 - _tmp268 * _tmp274;
  const Scalar _tmp276 = Scalar(0.71007031138673404) * _tmp272;
  const Scalar _tmp277 = _tmp275 * _tmp276;
  const Scalar _tmp278 = std::sinh(_tmp277);
  const Scalar _tmp279 = _tmp163 * _tmp81;
  const Scalar _tmp280 = _tmp132 * _tmp142 - _tmp133 * _tmp142 + _tmp154 * _tmp279 -
                         _tmp158 * _tmp163 - _tmp234 - _tmp238;
  const Scalar _tmp281 = std::pow(_tmp268, Scalar(-2));
  const Scalar _tmp282 = Scalar(0.71007031138673404) * _tmp281;
  const Scalar _tmp283 = _tmp280 * _tmp282;
  const Scalar _tmp284 = Scalar(1.4083112389913199) * _tmp268;
  const Scalar _tmp285 = _tmp271 * _tmp281;
  const Scalar _tmp286 =
      std::pow(Scalar(std::pow(_tmp271, Scalar(2)) * _tmp281 + 1), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp287 =
      _tmp286 *
      (_tmp272 * (-_tmp123 * _tmp270 - _tmp126 * _tmp270 + _tmp130 * _tmp270 - _tmp149 * _tmp269 -
                  _tmp151 * _tmp269 + _tmp153 * _tmp269 + _tmp170 * _tmp208 - _tmp170 * _tmp210) -
       _tmp280 * _tmp285);
  const Scalar _tmp288 = Scalar(1.0) / (_tmp267);
  const Scalar _tmp289 = Scalar(1.0) * _tmp273;
  const Scalar _tmp290 = Scalar(1.0) * std::sinh(_tmp289);
  const Scalar _tmp291 = Scalar(1.4083112389913199) * _tmp276 * p_a(2, 0) -
                         Scalar(1.4083112389913199) * std::cosh(_tmp277) +
                         Scalar(1.4083112389913199) * std::cosh(_tmp289);
  const Scalar _tmp292 = _tmp33 * _tmp86;
  const Scalar _tmp293 = -_tmp292 + _tmp35 + _tmp51 * _tmp88;
  const Scalar _tmp294 =
      _tmp107 * _tmp88 - _tmp108 * _tmp293 - _tmp17 * _tmp292 + _tmp17 * _tmp35 + _tmp30 * _tmp88;
  const Scalar _tmp295 = _tmp128 * _tmp294;
  const Scalar _tmp296 = -_tmp111 * _tmp88 + _tmp292 * _tmp41 - _tmp293 * _tmp70 - _tmp42;
  const Scalar _tmp297 = -_tmp293 * _tmp79 - _tmp296 * _tmp78 + _tmp67 * _tmp88 - _tmp91;
  const Scalar _tmp298 = _tmp297 * _tmp94;
  const Scalar _tmp299 = _tmp117 * (-_tmp146 * _tmp294 + _tmp150 * _tmp298);
  const Scalar _tmp300 = _tmp148 * _tmp297;
  const Scalar _tmp301 = _tmp152 * _tmp295;
  const Scalar _tmp302 = _tmp299 + _tmp300 - _tmp301;
  const Scalar _tmp303 = _tmp157 * _tmp298;
  const Scalar _tmp304 = -_tmp134 * _tmp302 - _tmp159 * _tmp293 + _tmp303 * _tmp54;
  const Scalar _tmp305 = _tmp119 * _tmp298;
  const Scalar _tmp306 = _tmp122 * _tmp297;
  const Scalar _tmp307 = _tmp117 * (-_tmp114 * _tmp294 + _tmp125 * _tmp298);
  const Scalar _tmp308 = _tmp127 * _tmp295;
  const Scalar _tmp309 = _tmp306 + _tmp307 - _tmp308;
  const Scalar _tmp310 = -_tmp120 * _tmp293 - _tmp134 * _tmp309 + _tmp305 * _tmp54;
  const Scalar _tmp311 =
      -_tmp105 * (_tmp101 * _tmp298 - _tmp293 * _tmp99 - _tmp298 * _tmp97 + _tmp84 * _tmp88) -
      _tmp143 * (_tmp120 * _tmp88 + _tmp155 * _tmp309 - _tmp305 * _tmp96 + _tmp310 * _tmp98) -
      _tmp164 * (_tmp155 * _tmp302 + _tmp159 * _tmp88 - _tmp303 * _tmp96 + _tmp304 * _tmp98) -
      _tmp171 * (_tmp165 * _tmp88 - _tmp166 * _tmp293 + _tmp167 * _tmp295 - _tmp168 * _tmp295);
  const Scalar _tmp312 =
      -_tmp165 * _tmp296 - _tmp203 * _tmp295 + _tmp204 * _tmp295 + _tmp205 * _tmp297;
  const Scalar _tmp313 = _tmp207 * _tmp295;
  const Scalar _tmp314 = _tmp209 * _tmp297;
  const Scalar _tmp315 = -_tmp159 * _tmp296 + _tmp180 * _tmp299 + _tmp180 * _tmp300 -
                         _tmp180 * _tmp301 - _tmp181 * _tmp302 + _tmp303 * _tmp73;
  const Scalar _tmp316 = _tmp296 * _tmp81;
  const Scalar _tmp317 = -_tmp119 * _tmp316 + _tmp180 * _tmp306 + _tmp180 * _tmp307 -
                         _tmp180 * _tmp308 - _tmp181 * _tmp309 + _tmp305 * _tmp73;
  const Scalar _tmp318 = _tmp212 * (_tmp211 * _tmp298 - _tmp316 * _tmp83);
  const Scalar _tmp319 =
      -_tmp200 * _tmp311 +
      _tmp201 * (_tmp184 * (_tmp176 * _tmp306 + _tmp176 * _tmp307 - _tmp176 * _tmp308 -
                            _tmp183 * _tmp317) +
                 _tmp188 * (_tmp176 * _tmp299 + _tmp176 * _tmp300 - _tmp176 * _tmp301 -
                            _tmp183 * _tmp315) +
                 _tmp194 * (-_tmp183 * _tmp312 - _tmp313 + _tmp314) - Scalar(1.0) * _tmp318);
  const Scalar _tmp320 = Scalar(1.4083112389913199) * _tmp311;
  const Scalar _tmp321 = _tmp233 * _tmp298;
  const Scalar _tmp322 = _tmp194 * _tmp295;
  const Scalar _tmp323 = -_tmp232 * _tmp293 + _tmp235 * _tmp321 + _tmp235 * _tmp322 -
                         _tmp237 * _tmp293 + _tmp239 * _tmp304 + _tmp240 * _tmp310;
  const Scalar _tmp324 = _tmp244 * _tmp323;
  const Scalar _tmp325 =
      _tmp255 * (_tmp251 * (_tmp247 * _tmp312 + _tmp248 * _tmp315 + _tmp249 * _tmp317 + _tmp318) -
                 _tmp256 * _tmp323);
  const Scalar _tmp326 = -_tmp142 * _tmp305 + _tmp142 * _tmp309 * _tmp81 - _tmp163 * _tmp303 +
                         _tmp279 * _tmp302 - _tmp321 - _tmp322;
  const Scalar _tmp327 = _tmp282 * _tmp326;
  const Scalar _tmp328 =
      _tmp286 *
      (_tmp272 * (_tmp170 * _tmp313 - _tmp170 * _tmp314 - _tmp269 * _tmp299 - _tmp269 * _tmp300 +
                  _tmp269 * _tmp301 - _tmp270 * _tmp306 - _tmp270 * _tmp307 + _tmp270 * _tmp308) -
       _tmp285 * _tmp326);

  // Output terms (1)
  Eigen::Matrix<Scalar, 4, 3> _res;

  _res(0, 0) = 0;
  _res(1, 0) =
      -_tmp222 *
          (-_tmp172 * _tmp221 + _tmp215 * _tmp219 -
           _tmp226 * (-_tmp172 * _tmp227 + _tmp224 * (-_tmp215 * _tmp228 - _tmp217 * _tmp229))) -
      _tmp229 * _tmp230;
  _res(2, 0) =
      -_tmp241 * _tmp264 -
      _tmp263 *
          (-_tmp245 * p_c(2, 0) + _tmp254 * _tmp257 -
           _tmp262 * (-_tmp245 * _tmp259 + _tmp260 * (-_tmp241 * _tmp258 - _tmp257 * _tmp263)));
  _res(3, 0) =
      -_tmp280 * _tmp291 -
      _tmp284 *
          (-_tmp278 * (-_tmp275 * _tmp283 +
                       _tmp276 * (-_tmp266 * _tmp288 - _tmp274 * _tmp280 - _tmp284 * _tmp287)) -
           _tmp283 * p_a(2, 0) + _tmp287 * _tmp290);
  _res(0, 1) = 0;
  _res(1, 1) =
      -_tmp222 *
          (_tmp219 * _tmp319 - _tmp221 * _tmp311 -
           _tmp226 * (_tmp224 * (-_tmp217 * _tmp320 - _tmp228 * _tmp319) - _tmp227 * _tmp311)) -
      _tmp230 * _tmp320;
  _res(2, 1) =
      -_tmp263 *
          (_tmp254 * _tmp325 -
           _tmp262 * (-_tmp259 * _tmp324 + _tmp260 * (-_tmp258 * _tmp323 - _tmp263 * _tmp325)) -
           _tmp324 * p_c(2, 0)) -
      _tmp264 * _tmp323;
  _res(3, 1) =
      -_tmp284 *
          (-_tmp278 * (-_tmp275 * _tmp327 +
                       _tmp276 * (-_tmp265 * _tmp288 - _tmp274 * _tmp326 - _tmp284 * _tmp328)) +
           _tmp290 * _tmp328 - _tmp327 * p_a(2, 0)) -
      _tmp291 * _tmp326;
  _res(0, 2) = 0;
  _res(1, 2) = 0;
  _res(2, 2) = 0;
  _res(3, 2) = Scalar(-1.0);

  return _res;
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym

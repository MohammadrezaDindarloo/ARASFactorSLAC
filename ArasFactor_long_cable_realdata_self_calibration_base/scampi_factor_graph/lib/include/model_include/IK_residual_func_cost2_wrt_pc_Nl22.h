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
 * Symbolic function: IK_residual_func_cost2_wrt_pc_Nl22
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
Eigen::Matrix<Scalar, 4, 3> IkResidualFuncCost2WrtPcNl22(
    const Scalar fh1, const Scalar fv1, const sym::Rot3<Scalar>& DeltaRot,
    const Eigen::Matrix<Scalar, 3, 1>& position_vector, const Eigen::Matrix<Scalar, 4, 1>& encoder,
    const Eigen::Matrix<Scalar, 3, 1>& p_a, const Eigen::Matrix<Scalar, 3, 1>& p_b,
    const Eigen::Matrix<Scalar, 3, 1>& p_c, const Eigen::Matrix<Scalar, 3, 1>& p_d,
    const sym::Rot3<Scalar>& Rot_init, const Scalar epsilon) {
  // Total ops: 1216

  // Unused inputs
  (void)encoder;
  (void)epsilon;

  // Input arrays
  const Eigen::Matrix<Scalar, 4, 1>& _DeltaRot = DeltaRot.Data();
  const Eigen::Matrix<Scalar, 4, 1>& _Rot_init = Rot_init.Data();

  // Intermediate terms (394)
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
  const Scalar _tmp7 = -2 * std::pow(_tmp0, Scalar(2));
  const Scalar _tmp8 = 1 - 2 * std::pow(_tmp3, Scalar(2));
  const Scalar _tmp9 = Scalar(0.20999999999999999) * _tmp7 + Scalar(0.20999999999999999) * _tmp8;
  const Scalar _tmp10 = 2 * _tmp3;
  const Scalar _tmp11 = _tmp1 * _tmp10;
  const Scalar _tmp12 = _tmp0 * _tmp4;
  const Scalar _tmp13 = _tmp11 + _tmp12;
  const Scalar _tmp14 = -Scalar(0.010999999999999999) * _tmp13;
  const Scalar _tmp15 = _tmp14 + _tmp9;
  const Scalar _tmp16 = _tmp15 + _tmp6;
  const Scalar _tmp17 = _tmp16 + position_vector(0, 0);
  const Scalar _tmp18 = _tmp17 - p_c(0, 0);
  const Scalar _tmp19 = Scalar(1.0) / (_tmp18);
  const Scalar _tmp20 = std::pow(_tmp18, Scalar(2));
  const Scalar _tmp21 = -2 * std::pow(_tmp1, Scalar(2));
  const Scalar _tmp22 = Scalar(0.20999999999999999) * _tmp21 + Scalar(0.20999999999999999) * _tmp8;
  const Scalar _tmp23 = _tmp0 * _tmp10;
  const Scalar _tmp24 = _tmp1 * _tmp4;
  const Scalar _tmp25 = _tmp23 - _tmp24;
  const Scalar _tmp26 = -Scalar(0.010999999999999999) * _tmp25;
  const Scalar _tmp27 = Scalar(0.20999999999999999) * _tmp2 + Scalar(0.20999999999999999) * _tmp5;
  const Scalar _tmp28 = _tmp26 + _tmp27;
  const Scalar _tmp29 = _tmp22 + _tmp28;
  const Scalar _tmp30 = _tmp29 + position_vector(1, 0);
  const Scalar _tmp31 = _tmp30 - p_c(1, 0);
  const Scalar _tmp32 = std::pow(_tmp31, Scalar(2));
  const Scalar _tmp33 = _tmp20 + _tmp32;
  const Scalar _tmp34 = std::sqrt(_tmp33);
  const Scalar _tmp35 = _tmp19 * _tmp34;
  const Scalar _tmp36 = Scalar(333.54000000000002) - fv1;
  const Scalar _tmp37 = -_tmp22;
  const Scalar _tmp38 = _tmp26 - _tmp27;
  const Scalar _tmp39 = _tmp37 + _tmp38;
  const Scalar _tmp40 = _tmp39 + position_vector(1, 0);
  const Scalar _tmp41 = _tmp40 - p_a(1, 0);
  const Scalar _tmp42 = -_tmp6;
  const Scalar _tmp43 = _tmp14 - _tmp9;
  const Scalar _tmp44 = _tmp42 + _tmp43;
  const Scalar _tmp45 = _tmp44 + position_vector(0, 0);
  const Scalar _tmp46 = _tmp45 - p_a(0, 0);
  const Scalar _tmp47 = std::pow(_tmp46, Scalar(2));
  const Scalar _tmp48 = std::pow(_tmp41, Scalar(2)) + _tmp47;
  const Scalar _tmp49 = std::pow(_tmp48, Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp50 = _tmp41 * _tmp49;
  const Scalar _tmp51 = _tmp19 * _tmp31;
  const Scalar _tmp52 = _tmp46 * _tmp49;
  const Scalar _tmp53 = _tmp51 * _tmp52;
  const Scalar _tmp54 = -_tmp50 + _tmp53;
  const Scalar _tmp55 = Scalar(1.0) / (_tmp54);
  const Scalar _tmp56 = _tmp52 * _tmp55;
  const Scalar _tmp57 = Scalar(1.0) * _tmp29;
  const Scalar _tmp58 = -_tmp57;
  const Scalar _tmp59 = Scalar(1.0) / (_tmp39 + _tmp58);
  const Scalar _tmp60 = Scalar(1.0) * _tmp16;
  const Scalar _tmp61 = -_tmp44 + _tmp60;
  const Scalar _tmp62 = _tmp59 * _tmp61;
  const Scalar _tmp63 = _tmp57 * _tmp62 + _tmp60;
  const Scalar _tmp64 = 0;
  const Scalar _tmp65 = Scalar(0.20999999999999999) * _tmp11 - Scalar(0.20999999999999999) * _tmp12;
  const Scalar _tmp66 = -Scalar(0.010999999999999999) * _tmp21 -
                        Scalar(0.010999999999999999) * _tmp7 + Scalar(-0.010999999999999999);
  const Scalar _tmp67 = Scalar(0.20999999999999999) * _tmp23 + Scalar(0.20999999999999999) * _tmp24;
  const Scalar _tmp68 = _tmp66 - _tmp67;
  const Scalar _tmp69 = _tmp65 + _tmp68;
  const Scalar _tmp70 = _tmp15 + _tmp42;
  const Scalar _tmp71 = _tmp70 + position_vector(0, 0);
  const Scalar _tmp72 = _tmp71 - p_b(0, 0);
  const Scalar _tmp73 = _tmp28 + _tmp37;
  const Scalar _tmp74 = _tmp73 + position_vector(1, 0);
  const Scalar _tmp75 = _tmp74 - p_b(1, 0);
  const Scalar _tmp76 = std::pow(Scalar(std::pow(_tmp72, Scalar(2)) + std::pow(_tmp75, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp77 = _tmp72 * _tmp76;
  const Scalar _tmp78 = -_tmp65;
  const Scalar _tmp79 = _tmp68 + _tmp78;
  const Scalar _tmp80 = _tmp66 + _tmp67;
  const Scalar _tmp81 = _tmp65 + _tmp80;
  const Scalar _tmp82 = _tmp52 * _tmp81;
  const Scalar _tmp83 = -_tmp52 * _tmp79 + _tmp82;
  const Scalar _tmp84 = _tmp75 * _tmp76;
  const Scalar _tmp85 = _tmp51 * _tmp77 - _tmp84;
  const Scalar _tmp86 = _tmp55 * _tmp85;
  const Scalar _tmp87 = _tmp77 * _tmp81;
  const Scalar _tmp88 = _tmp51 * _tmp81;
  const Scalar _tmp89 = _tmp50 * _tmp79 - _tmp52 * _tmp88;
  const Scalar _tmp90 = _tmp69 * _tmp84 - _tmp77 * _tmp88 - _tmp86 * _tmp89;
  const Scalar _tmp91 = -_tmp62 * _tmp90 - _tmp69 * _tmp77 - _tmp83 * _tmp86 + _tmp87;
  const Scalar _tmp92 = Scalar(1.0) / (_tmp91);
  const Scalar _tmp93 = _tmp85 * _tmp92;
  const Scalar _tmp94 = _tmp64 * _tmp93;
  const Scalar _tmp95 = _tmp77 * _tmp92;
  const Scalar _tmp96 = _tmp64 * _tmp95;
  const Scalar _tmp97 = _tmp36 * (-_tmp56 * _tmp94 + _tmp96);
  const Scalar _tmp98 = Scalar(1.0) / (_tmp34);
  const Scalar _tmp99 = _tmp29 * _tmp98;
  const Scalar _tmp100 = _tmp16 * _tmp98;
  const Scalar _tmp101 = _tmp100 * _tmp31 - _tmp18 * _tmp99;
  const Scalar _tmp102 = _tmp101 * _tmp35;
  const Scalar _tmp103 = _tmp102 * _tmp52 + _tmp39 * _tmp52 - _tmp44 * _tmp50;
  const Scalar _tmp104 = Scalar(1.0) * _tmp55;
  const Scalar _tmp105 = Scalar(1.0) * _tmp59;
  const Scalar _tmp106 = _tmp105 * _tmp61;
  const Scalar _tmp107 = _tmp106 * _tmp89;
  const Scalar _tmp108 = -_tmp104 * _tmp83 + _tmp107 * _tmp55;
  const Scalar _tmp109 = _tmp102 * _tmp77 - _tmp103 * _tmp86 - _tmp70 * _tmp84 + _tmp73 * _tmp77;
  const Scalar _tmp110 = _tmp109 * _tmp92;
  const Scalar _tmp111 = -_tmp103 * _tmp104 - _tmp108 * _tmp110;
  const Scalar _tmp112 = Scalar(1.0) / (_tmp109);
  const Scalar _tmp113 = _tmp111 * _tmp112;
  const Scalar _tmp114 = _tmp113 * _tmp91;
  const Scalar _tmp115 = _tmp108 + _tmp114;
  const Scalar _tmp116 = _tmp115 * _tmp95;
  const Scalar _tmp117 = -_tmp115 * _tmp93 + Scalar(1.0);
  const Scalar _tmp118 = _tmp43 + _tmp6;
  const Scalar _tmp119 = _tmp118 - p_d(0, 0) + position_vector(0, 0);
  const Scalar _tmp120 = _tmp22 + _tmp38;
  const Scalar _tmp121 = _tmp120 - p_d(1, 0) + position_vector(1, 0);
  const Scalar _tmp122 =
      std::pow(Scalar(std::pow(_tmp119, Scalar(2)) + std::pow(_tmp121, Scalar(2))),
               Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp123 = _tmp121 * _tmp122;
  const Scalar _tmp124 = _tmp123 * fh1;
  const Scalar _tmp125 = _tmp124 * (_tmp116 + _tmp117 * _tmp56);
  const Scalar _tmp126 = Scalar(1.0) * _tmp112;
  const Scalar _tmp127 = _tmp52 * _tmp86;
  const Scalar _tmp128 = _tmp126 * _tmp77;
  const Scalar _tmp129 = _tmp119 * _tmp122;
  const Scalar _tmp130 = fh1 * (_tmp118 * _tmp123 - _tmp120 * _tmp129);
  const Scalar _tmp131 = _tmp130 * (-_tmp126 * _tmp127 + _tmp128);
  const Scalar _tmp132 = _tmp51 * _tmp55;
  const Scalar _tmp133 = _tmp132 * _tmp89 + _tmp88;
  const Scalar _tmp134 = _tmp132 * _tmp83 - _tmp133 * _tmp62 - _tmp81;
  const Scalar _tmp135 = _tmp103 * _tmp55;
  const Scalar _tmp136 = -_tmp102 - _tmp110 * _tmp134 + _tmp135 * _tmp51;
  const Scalar _tmp137 = _tmp112 * _tmp91;
  const Scalar _tmp138 = _tmp136 * _tmp137;
  const Scalar _tmp139 = _tmp134 + _tmp138;
  const Scalar _tmp140 = _tmp139 * _tmp95;
  const Scalar _tmp141 = -_tmp139 * _tmp93 - _tmp51;
  const Scalar _tmp142 = _tmp129 * fh1;
  const Scalar _tmp143 = _tmp142 * (_tmp140 + _tmp141 * _tmp56 + Scalar(1.0));
  const Scalar _tmp144 = -_tmp125 * _tmp35 - _tmp131 * _tmp35 - _tmp143 * _tmp35 - _tmp35 * _tmp97;
  const Scalar _tmp145 = Scalar(1.0) / (_tmp144);
  const Scalar _tmp146 = std::pow(_tmp54, Scalar(-2));
  const Scalar _tmp147 = Scalar(1.0) / (_tmp20);
  const Scalar _tmp148 = _tmp147 * _tmp31;
  const Scalar _tmp149 = _tmp146 * _tmp148;
  const Scalar _tmp150 = _tmp52 * _tmp83;
  const Scalar _tmp151 = _tmp149 * _tmp150;
  const Scalar _tmp152 = _tmp52 * _tmp89;
  const Scalar _tmp153 = _tmp149 * _tmp152;
  const Scalar _tmp154 = _tmp148 * _tmp81;
  const Scalar _tmp155 = _tmp148 * _tmp55;
  const Scalar _tmp156 = _tmp155 * _tmp89;
  const Scalar _tmp157 = _tmp127 * _tmp154 - _tmp148 * _tmp87 + _tmp153 * _tmp85 - _tmp156 * _tmp77;
  const Scalar _tmp158 = _tmp155 * _tmp83;
  const Scalar _tmp159 = _tmp151 * _tmp85 - _tmp157 * _tmp62 - _tmp158 * _tmp77;
  const Scalar _tmp160 = _tmp58 + _tmp73;
  const Scalar _tmp161 = _tmp160 * _tmp62;
  const Scalar _tmp162 = Scalar(1.0) / (-_tmp161 + _tmp60 - _tmp70);
  const Scalar _tmp163 = Scalar(1.0) * _tmp162;
  const Scalar _tmp164 = _tmp112 * _tmp163;
  const Scalar _tmp165 = _tmp159 * _tmp164;
  const Scalar _tmp166 = _tmp160 * _tmp163;
  const Scalar _tmp167 = _tmp135 * _tmp148;
  const Scalar _tmp168 = std::pow(_tmp33, Scalar(Scalar(-3) / Scalar(2)));
  const Scalar _tmp169 = _tmp168 * _tmp29;
  const Scalar _tmp170 = _tmp16 * _tmp168;
  const Scalar _tmp171 = _tmp18 * _tmp31;
  const Scalar _tmp172 = _tmp35 * (-_tmp169 * _tmp20 + _tmp170 * _tmp171 + _tmp99);
  const Scalar _tmp173 = _tmp147 * _tmp34;
  const Scalar _tmp174 = _tmp101 * _tmp173;
  const Scalar _tmp175 = _tmp103 * _tmp52;
  const Scalar _tmp176 = _tmp149 * _tmp175;
  const Scalar _tmp177 = _tmp101 * _tmp98;
  const Scalar _tmp178 = _tmp172 * _tmp52 + _tmp174 * _tmp52 - _tmp177 * _tmp52;
  const Scalar _tmp179 = _tmp177 * _tmp77;
  const Scalar _tmp180 = -_tmp167 * _tmp77 + _tmp172 * _tmp77 + _tmp174 * _tmp77 +
                         _tmp176 * _tmp85 - _tmp178 * _tmp86 - _tmp179;
  const Scalar _tmp181 = std::pow(_tmp109, Scalar(-2));
  const Scalar _tmp182 = _tmp180 * _tmp181;
  const Scalar _tmp183 = _tmp182 * _tmp91;
  const Scalar _tmp184 = Scalar(1.0) * _tmp90;
  const Scalar _tmp185 = _tmp112 * _tmp166;
  const Scalar _tmp186 =
      -_tmp126 * _tmp157 + _tmp159 * _tmp185 - _tmp166 * _tmp183 + _tmp182 * _tmp184;
  const Scalar _tmp187 = _tmp163 * _tmp183;
  const Scalar _tmp188 = Scalar(1.0) * _tmp130;
  const Scalar _tmp189 = std::pow(_tmp91, Scalar(-2));
  const Scalar _tmp190 = _tmp189 * _tmp90;
  const Scalar _tmp191 = _tmp159 * _tmp190;
  const Scalar _tmp192 = _tmp157 * _tmp92;
  const Scalar _tmp193 = _tmp191 * _tmp64 - _tmp192 * _tmp64;
  const Scalar _tmp194 = _tmp105 * _tmp36;
  const Scalar _tmp195 = _tmp112 * _tmp136;
  const Scalar _tmp196 = _tmp159 * _tmp195;
  const Scalar _tmp197 = _tmp32 / [&]() {
    const Scalar base = _tmp18;
    return base * base * base;
  }();
  const Scalar _tmp198 = _tmp146 * _tmp197;
  const Scalar _tmp199 = -_tmp152 * _tmp198 + _tmp154 + _tmp156 - _tmp197 * _tmp55 * _tmp82;
  const Scalar _tmp200 = -_tmp150 * _tmp198 + _tmp158 - _tmp199 * _tmp62;
  const Scalar _tmp201 = _tmp159 * _tmp189;
  const Scalar _tmp202 = _tmp109 * _tmp134;
  const Scalar _tmp203 = _tmp134 * _tmp92;
  const Scalar _tmp204 =
      _tmp137 * (-_tmp110 * _tmp200 + _tmp132 * _tmp178 + _tmp167 - _tmp172 - _tmp174 -
                 _tmp175 * _tmp198 + _tmp177 - _tmp180 * _tmp203 + _tmp201 * _tmp202);
  const Scalar _tmp205 = _tmp136 * _tmp183;
  const Scalar _tmp206 = _tmp196 + _tmp200 + _tmp204 - _tmp205;
  const Scalar _tmp207 = _tmp90 * _tmp92;
  const Scalar _tmp208 = _tmp160 * _tmp162;
  const Scalar _tmp209 = _tmp139 * _tmp191 - _tmp139 * _tmp192 + _tmp196 * _tmp208 + _tmp199 +
                         _tmp204 * _tmp208 - _tmp205 * _tmp208 - _tmp206 * _tmp207;
  const Scalar _tmp210 = Scalar(1.0) * _tmp142;
  const Scalar _tmp211 = _tmp115 * _tmp159;
  const Scalar _tmp212 = _tmp104 * _tmp52;
  const Scalar _tmp213 = _tmp111 * _tmp183;
  const Scalar _tmp214 = _tmp113 * _tmp159;
  const Scalar _tmp215 = _tmp154 * _tmp56;
  const Scalar _tmp216 = _tmp149 * _tmp52;
  const Scalar _tmp217 = -_tmp106 * _tmp215 - _tmp107 * _tmp216 + Scalar(1.0) * _tmp151;
  const Scalar _tmp218 = _tmp108 * _tmp92;
  const Scalar _tmp219 = _tmp108 * _tmp109;
  const Scalar _tmp220 = _tmp137 * (-_tmp104 * _tmp178 - _tmp110 * _tmp217 + Scalar(1.0) * _tmp176 -
                                    _tmp180 * _tmp218 + _tmp201 * _tmp219);
  const Scalar _tmp221 = -_tmp213 + _tmp214 + _tmp217 + _tmp220;
  const Scalar _tmp222 = -_tmp115 * _tmp192 + Scalar(1.0) * _tmp153 + _tmp154 * _tmp212 +
                         _tmp190 * _tmp211 - _tmp207 * _tmp221 - _tmp208 * _tmp213 +
                         _tmp208 * _tmp214 + _tmp208 * _tmp220;
  const Scalar _tmp223 = Scalar(1.0) * _tmp124;
  const Scalar _tmp224 = _tmp125 * _tmp98;
  const Scalar _tmp225 = _tmp47 / _tmp48;
  const Scalar _tmp226 = _tmp149 * _tmp225;
  const Scalar _tmp227 = _tmp127 * _tmp64;
  const Scalar _tmp228 = _tmp64 * _tmp77;
  const Scalar _tmp229 = _tmp35 * _tmp36;
  const Scalar _tmp230 = _tmp112 * _tmp148 * _tmp77;
  const Scalar _tmp231 = Scalar(1.0) * _tmp127;
  const Scalar _tmp232 = Scalar(1.0) * _tmp77;
  const Scalar _tmp233 = _tmp130 * _tmp35;
  const Scalar _tmp234 = _tmp143 * _tmp98;
  const Scalar _tmp235 = _tmp189 * _tmp85;
  const Scalar _tmp236 = -_tmp116 * _tmp148 + _tmp211 * _tmp235 - _tmp221 * _tmp93;
  const Scalar _tmp237 = _tmp189 * _tmp77;
  const Scalar _tmp238 = _tmp124 * _tmp35;
  const Scalar _tmp239 = _tmp141 * _tmp225;
  const Scalar _tmp240 = _tmp139 * _tmp235;
  const Scalar _tmp241 = -_tmp140 * _tmp148 - _tmp148 + _tmp159 * _tmp240 - _tmp206 * _tmp93;
  const Scalar _tmp242 = _tmp139 * _tmp77;
  const Scalar _tmp243 = _tmp142 * _tmp35;
  const Scalar _tmp244 = _tmp131 * _tmp98;
  const Scalar _tmp245 = _tmp97 * _tmp98;
  const Scalar _tmp246 =
      -_tmp125 * _tmp173 - _tmp131 * _tmp173 - _tmp143 * _tmp173 - _tmp173 * _tmp97 + _tmp224 -
      _tmp229 *
          (-_tmp148 * _tmp56 * _tmp96 + _tmp201 * _tmp227 - _tmp201 * _tmp228 + _tmp226 * _tmp94) -
      _tmp233 *
          (_tmp126 * _tmp226 * _tmp85 + _tmp182 * _tmp231 - _tmp182 * _tmp232 - _tmp212 * _tmp230) +
      _tmp234 -
      _tmp238 * (-_tmp117 * _tmp226 - _tmp211 * _tmp237 + _tmp221 * _tmp95 + _tmp236 * _tmp56) -
      _tmp243 * (-_tmp149 * _tmp239 - _tmp201 * _tmp242 + _tmp206 * _tmp95 + _tmp241 * _tmp56) +
      _tmp244 + _tmp245;
  const Scalar _tmp247 = -_tmp104 * _tmp89 + _tmp114 * _tmp208 - _tmp115 * _tmp207;
  const Scalar _tmp248 = fh1 * (_tmp78 + _tmp80);
  const Scalar _tmp249 = -_tmp120 * fv1 - _tmp123 * _tmp248 - Scalar(40.024799999999999) * _tmp25;
  const Scalar _tmp250 = _tmp161 * _tmp163 + Scalar(1.0);
  const Scalar _tmp251 = _tmp163 * _tmp62;
  const Scalar _tmp252 = _tmp133 + _tmp138 * _tmp208 - _tmp139 * _tmp207;
  const Scalar _tmp253 = _tmp137 * _tmp163;
  const Scalar _tmp254 = -_tmp126 * _tmp90 + _tmp137 * _tmp166;
  const Scalar _tmp255 = _tmp162 * _tmp63;
  const Scalar _tmp256 = -_tmp160 * _tmp255 - _tmp207 * _tmp64 + _tmp58;
  const Scalar _tmp257 = _tmp118 * fv1 + _tmp129 * _tmp248 + Scalar(40.024799999999999) * _tmp13;
  const Scalar _tmp258 =
      _tmp188 * (-_tmp105 * _tmp254 + _tmp253) +
      _tmp210 * (-_tmp105 * _tmp252 + _tmp138 * _tmp163) +
      _tmp223 * (-_tmp105 * _tmp247 + _tmp114 * _tmp163) +
      Scalar(1.0) * _tmp249 * (-_tmp105 * _tmp250 + _tmp251) +
      Scalar(1.0) * _tmp257 * (-_tmp163 + _tmp166 * _tmp59) +
      Scalar(1.0) * _tmp36 * (-_tmp105 * _tmp256 - _tmp163 * _tmp63 + Scalar(1.0));
  const Scalar _tmp259 = std::pow(_tmp144, Scalar(-2));
  const Scalar _tmp260 = _tmp258 * _tmp259;
  const Scalar _tmp261 =
      _tmp145 * (_tmp188 * (-_tmp105 * _tmp186 + _tmp165 - _tmp187) - _tmp193 * _tmp194 +
                 _tmp210 * (-_tmp105 * _tmp209 + _tmp163 * _tmp196 + _tmp163 * _tmp204 -
                            _tmp163 * _tmp205) +
                 _tmp223 * (-_tmp105 * _tmp222 - _tmp163 * _tmp213 + _tmp163 * _tmp214 +
                            _tmp163 * _tmp220)) -
      _tmp246 * _tmp260;
  const Scalar _tmp262 =
      std::pow(Scalar(std::pow(_tmp258, Scalar(2)) * _tmp259 + 1), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp263 = std::asinh(_tmp145 * _tmp258);
  const Scalar _tmp264 = Scalar(1.0) * _tmp263;
  const Scalar _tmp265 = Scalar(1.0) * _tmp262 * std::cosh(_tmp264);
  const Scalar _tmp266 = Scalar(1.4083112389913199) * _tmp144;
  const Scalar _tmp267 = -_tmp30 + p_c(1, 0);
  const Scalar _tmp268 = -_tmp17 + p_c(0, 0);
  const Scalar _tmp269 = std::pow(_tmp267, Scalar(2)) + std::pow(_tmp268, Scalar(2));
  const Scalar _tmp270 = std::sqrt(_tmp269);
  const Scalar _tmp271 = -_tmp263 * _tmp266 - _tmp270;
  const Scalar _tmp272 = Scalar(0.71007031138673404) * _tmp145;
  const Scalar _tmp273 = _tmp271 * _tmp272;
  const Scalar _tmp274 = std::cosh(_tmp273);
  const Scalar _tmp275 = Scalar(1.0) / (_tmp270);
  const Scalar _tmp276 = Scalar(1.4083112389913199) * _tmp263;
  const Scalar _tmp277 = _tmp262 * _tmp266;
  const Scalar _tmp278 = Scalar(0.71007031138673404) * _tmp259 * _tmp271;
  const Scalar _tmp279 = -_tmp81 + p_c(2, 0) - position_vector(2, 0);
  const Scalar _tmp280 =
      std::pow(Scalar(_tmp269 + std::pow(_tmp279, Scalar(2))), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp281 = -Scalar(1.4083112389913199) * std::sinh(_tmp264) -
                         Scalar(1.4083112389913199) * std::sinh(_tmp273);
  const Scalar _tmp282 = _tmp124 * _tmp59;
  const Scalar _tmp283 = _tmp142 * _tmp59;
  const Scalar _tmp284 = _tmp36 * _tmp59;
  const Scalar _tmp285 = _tmp130 * _tmp59;
  const Scalar _tmp286 = _tmp163 * _tmp257;
  const Scalar _tmp287 = -_tmp160 * _tmp286 * _tmp59 + _tmp247 * _tmp282 +
                         _tmp249 * _tmp250 * _tmp59 + _tmp252 * _tmp283 + _tmp254 * _tmp285 +
                         _tmp256 * _tmp284;
  const Scalar _tmp288 = _tmp55 * fh1;
  const Scalar _tmp289 = _tmp129 * _tmp288;
  const Scalar _tmp290 = _tmp126 * _tmp130;
  const Scalar _tmp291 = _tmp123 * _tmp288;
  const Scalar _tmp292 = _tmp36 * _tmp64;
  const Scalar _tmp293 = _tmp292 * _tmp93;
  const Scalar _tmp294 =
      _tmp117 * _tmp291 + _tmp141 * _tmp289 - _tmp290 * _tmp86 - _tmp293 * _tmp55;
  const Scalar _tmp295 = Scalar(1.0) / (_tmp294);
  const Scalar _tmp296 = std::asinh(_tmp287 * _tmp295);
  const Scalar _tmp297 = Scalar(1.0) * _tmp296;
  const Scalar _tmp298 = Scalar(1.0) * std::cosh(_tmp297);
  const Scalar _tmp299 = std::pow(_tmp294, Scalar(-2));
  const Scalar _tmp300 =
      std::pow(Scalar(std::pow(_tmp287, Scalar(2)) * _tmp299 + 1), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp301 = _tmp182 * _tmp188;
  const Scalar _tmp302 = _tmp141 * _tmp142;
  const Scalar _tmp303 = _tmp292 * _tmp92;
  const Scalar _tmp304 = _tmp303 * _tmp77;
  const Scalar _tmp305 = _tmp124 * _tmp52;
  const Scalar _tmp306 = _tmp201 * _tmp292;
  const Scalar _tmp307 = -_tmp104 * _tmp130 * _tmp230 - _tmp117 * _tmp149 * _tmp305 -
                         _tmp155 * _tmp304 + _tmp216 * _tmp290 * _tmp85 + _tmp216 * _tmp293 -
                         _tmp216 * _tmp302 + _tmp236 * _tmp291 + _tmp241 * _tmp289 +
                         _tmp301 * _tmp86 + _tmp306 * _tmp86;
  const Scalar _tmp308 = _tmp287 * _tmp299;
  const Scalar _tmp309 =
      _tmp300 *
      (_tmp295 * (_tmp186 * _tmp285 + _tmp193 * _tmp284 + _tmp209 * _tmp283 + _tmp222 * _tmp282) -
       _tmp307 * _tmp308);
  const Scalar _tmp310 = Scalar(1.4083112389913199) * _tmp294;
  const Scalar _tmp311 = Scalar(1.4083112389913199) * _tmp307;
  const Scalar _tmp312 = Scalar(0.71007031138673404) * _tmp295;
  const Scalar _tmp313 =
      -_tmp296 * _tmp310 - std::sqrt(Scalar(std::pow(Scalar(-_tmp40 + p_a(1, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp45 + p_a(0, 0)), Scalar(2))));
  const Scalar _tmp314 = Scalar(0.71007031138673404) * _tmp299 * _tmp313;
  const Scalar _tmp315 = _tmp312 * _tmp313;
  const Scalar _tmp316 = std::cosh(_tmp315);
  const Scalar _tmp317 = -std::sinh(_tmp297) - std::sinh(_tmp315);
  const Scalar _tmp318 = _tmp142 * _tmp92;
  const Scalar _tmp319 = _tmp124 * _tmp189;
  const Scalar _tmp320 = _tmp124 * _tmp92;
  const Scalar _tmp321 = _tmp139 * _tmp142;
  const Scalar _tmp322 = -_tmp201 * _tmp321 + _tmp206 * _tmp318 - _tmp211 * _tmp319 +
                         _tmp221 * _tmp320 - _tmp301 - _tmp306;
  const Scalar _tmp323 = _tmp115 * _tmp320 + _tmp139 * _tmp318 + _tmp290 + _tmp303;
  const Scalar _tmp324 = Scalar(1.0) / (_tmp323);
  const Scalar _tmp325 = _tmp124 * _tmp162;
  const Scalar _tmp326 = _tmp142 * _tmp162;
  const Scalar _tmp327 = -_tmp114 * _tmp325 - _tmp130 * _tmp253 - _tmp138 * _tmp326 -
                         _tmp249 * _tmp251 + _tmp255 * _tmp36 + _tmp286;
  const Scalar _tmp328 = std::asinh(_tmp324 * _tmp327);
  const Scalar _tmp329 = Scalar(1.4083112389913199) * _tmp328;
  const Scalar _tmp330 = Scalar(1.4083112389913199) * _tmp323;
  const Scalar _tmp331 = std::pow(_tmp323, Scalar(-2));
  const Scalar _tmp332 = _tmp327 * _tmp331;
  const Scalar _tmp333 =
      std::pow(Scalar(std::pow(_tmp327, Scalar(2)) * _tmp331 + 1), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp334 =
      _tmp333 *
      (-_tmp322 * _tmp332 +
       _tmp324 * (-_tmp130 * _tmp165 + _tmp130 * _tmp187 - _tmp196 * _tmp326 - _tmp204 * _tmp326 +
                  _tmp205 * _tmp326 + _tmp213 * _tmp325 - _tmp214 * _tmp325 - _tmp220 * _tmp325));
  const Scalar _tmp335 = Scalar(0.71007031138673404) * _tmp324;
  const Scalar _tmp336 =
      -_tmp328 * _tmp330 - std::sqrt(Scalar(std::pow(Scalar(-_tmp71 + p_b(0, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp74 + p_b(1, 0)), Scalar(2))));
  const Scalar _tmp337 = Scalar(0.71007031138673404) * _tmp331 * _tmp336;
  const Scalar _tmp338 = _tmp335 * _tmp336;
  const Scalar _tmp339 = std::cosh(_tmp338);
  const Scalar _tmp340 = Scalar(1.0) * _tmp328;
  const Scalar _tmp341 = Scalar(1.0) * std::cosh(_tmp340);
  const Scalar _tmp342 = -Scalar(1.4083112389913199) * std::sinh(_tmp338) -
                         Scalar(1.4083112389913199) * std::sinh(_tmp340);
  const Scalar _tmp343 = _tmp146 * _tmp19;
  const Scalar _tmp344 = _tmp19 * _tmp55;
  const Scalar _tmp345 = _tmp19 * _tmp81;
  const Scalar _tmp346 = _tmp345 * _tmp52;
  const Scalar _tmp347 = _tmp344 * _tmp89;
  const Scalar _tmp348 = _tmp343 * _tmp85;
  const Scalar _tmp349 =
      -_tmp152 * _tmp348 + _tmp345 * _tmp77 - _tmp346 * _tmp86 + _tmp347 * _tmp77;
  const Scalar _tmp350 = _tmp344 * _tmp83;
  const Scalar _tmp351 = -_tmp150 * _tmp348 - _tmp349 * _tmp62 + _tmp350 * _tmp77;
  const Scalar _tmp352 = _tmp189 * _tmp351;
  const Scalar _tmp353 = _tmp35 * (-_tmp100 - _tmp169 * _tmp171 + _tmp170 * _tmp32);
  const Scalar _tmp354 = -_tmp177 * _tmp53 + _tmp353 * _tmp52;
  const Scalar _tmp355 = _tmp103 * _tmp344;
  const Scalar _tmp356 = -_tmp175 * _tmp348 - _tmp179 * _tmp51 + _tmp353 * _tmp77 -
                         _tmp354 * _tmp86 + _tmp355 * _tmp77;
  const Scalar _tmp357 = _tmp181 * _tmp356;
  const Scalar _tmp358 = _tmp128 * _tmp344;
  const Scalar _tmp359 = _tmp357 * _tmp91;
  const Scalar _tmp360 = _tmp136 * _tmp359;
  const Scalar _tmp361 = _tmp195 * _tmp351;
  const Scalar _tmp362 = _tmp153 + _tmp215 - _tmp345 - _tmp347;
  const Scalar _tmp363 = _tmp151 - _tmp350 - _tmp362 * _tmp62;
  const Scalar _tmp364 =
      _tmp137 * (-_tmp110 * _tmp363 + _tmp132 * _tmp354 + _tmp176 + _tmp177 * _tmp51 +
                 _tmp202 * _tmp352 - _tmp203 * _tmp356 - _tmp353 - _tmp355);
  const Scalar _tmp365 = -_tmp360 + _tmp361 + _tmp363 + _tmp364;
  const Scalar _tmp366 = _tmp140 * _tmp19 + _tmp19 + _tmp240 * _tmp351 - _tmp365 * _tmp93;
  const Scalar _tmp367 = Scalar(1.0) * _tmp343;
  const Scalar _tmp368 = _tmp343 * _tmp52;
  const Scalar _tmp369 = _tmp106 * _tmp345 * _tmp56 + _tmp107 * _tmp368 - _tmp150 * _tmp367;
  const Scalar _tmp370 = _tmp137 * (-_tmp104 * _tmp354 - _tmp110 * _tmp369 - _tmp175 * _tmp367 -
                                    _tmp218 * _tmp356 + _tmp219 * _tmp352);
  const Scalar _tmp371 = _tmp113 * _tmp351;
  const Scalar _tmp372 = _tmp111 * _tmp359;
  const Scalar _tmp373 = _tmp369 + _tmp370 + _tmp371 - _tmp372;
  const Scalar _tmp374 = _tmp117 * _tmp343;
  const Scalar _tmp375 = _tmp115 * _tmp351;
  const Scalar _tmp376 = _tmp116 * _tmp19 + _tmp235 * _tmp375 - _tmp373 * _tmp93;
  const Scalar _tmp377 =
      _tmp224 * _tmp51 -
      _tmp229 * (-_tmp225 * _tmp343 * _tmp94 + _tmp227 * _tmp352 - _tmp228 * _tmp352 +
                 _tmp344 * _tmp52 * _tmp96) -
      _tmp233 * (-_tmp126 * _tmp225 * _tmp348 + _tmp231 * _tmp357 - _tmp232 * _tmp357 +
                 _tmp358 * _tmp52) +
      _tmp234 * _tmp51 -
      _tmp238 * (_tmp225 * _tmp374 - _tmp237 * _tmp375 + _tmp373 * _tmp95 + _tmp376 * _tmp56) -
      _tmp243 * (_tmp239 * _tmp343 - _tmp242 * _tmp352 + _tmp365 * _tmp95 + _tmp366 * _tmp56) +
      _tmp244 * _tmp51 + _tmp245 * _tmp51;
  const Scalar _tmp378 = _tmp349 * _tmp92;
  const Scalar _tmp379 = _tmp352 * _tmp90;
  const Scalar _tmp380 = -_tmp378 * _tmp64 + _tmp379 * _tmp64;
  const Scalar _tmp381 = _tmp163 * _tmp359;
  const Scalar _tmp382 = _tmp164 * _tmp351;
  const Scalar _tmp383 =
      -_tmp126 * _tmp349 - _tmp166 * _tmp359 + _tmp184 * _tmp357 + _tmp185 * _tmp351;
  const Scalar _tmp384 = -_tmp104 * _tmp346 - _tmp115 * _tmp378 - _tmp152 * _tmp367 +
                         _tmp190 * _tmp375 - _tmp207 * _tmp373 + _tmp208 * _tmp370 +
                         _tmp208 * _tmp371 - _tmp208 * _tmp372;
  const Scalar _tmp385 = -_tmp139 * _tmp378 + _tmp139 * _tmp379 - _tmp207 * _tmp365 -
                         _tmp208 * _tmp360 + _tmp208 * _tmp361 + _tmp208 * _tmp364 + _tmp362;
  const Scalar _tmp386 =
      _tmp145 * (_tmp188 * (-_tmp105 * _tmp383 - _tmp381 + _tmp382) - _tmp194 * _tmp380 +
                 _tmp210 * (-_tmp105 * _tmp385 - _tmp163 * _tmp360 + _tmp163 * _tmp361 +
                            _tmp163 * _tmp364) +
                 _tmp223 * (-_tmp105 * _tmp384 + _tmp163 * _tmp370 + _tmp163 * _tmp371 -
                            _tmp163 * _tmp372)) -
      _tmp260 * _tmp377;
  const Scalar _tmp387 = _tmp188 * _tmp357;
  const Scalar _tmp388 = _tmp292 * _tmp352;
  const Scalar _tmp389 = _tmp130 * _tmp358 + _tmp289 * _tmp366 - _tmp290 * _tmp348 * _tmp52 +
                         _tmp291 * _tmp376 - _tmp293 * _tmp368 + _tmp302 * _tmp368 +
                         _tmp304 * _tmp344 + _tmp305 * _tmp374 + _tmp387 * _tmp86 +
                         _tmp388 * _tmp86;
  const Scalar _tmp390 = Scalar(1.4083112389913199) * _tmp389;
  const Scalar _tmp391 =
      _tmp300 *
      (_tmp295 * (_tmp282 * _tmp384 + _tmp283 * _tmp385 + _tmp284 * _tmp380 + _tmp285 * _tmp383) -
       _tmp308 * _tmp389);
  const Scalar _tmp392 = _tmp318 * _tmp365 - _tmp319 * _tmp375 + _tmp320 * _tmp373 -
                         _tmp321 * _tmp352 - _tmp387 - _tmp388;
  const Scalar _tmp393 =
      _tmp333 *
      (_tmp324 * (_tmp130 * _tmp381 - _tmp130 * _tmp382 - _tmp325 * _tmp370 - _tmp325 * _tmp371 +
                  _tmp325 * _tmp372 + _tmp326 * _tmp360 - _tmp326 * _tmp361 - _tmp326 * _tmp364) -
       _tmp332 * _tmp392);

  // Output terms (1)
  Eigen::Matrix<Scalar, 4, 3> _res;

  _res(0, 0) = 0;
  _res(1, 0) =
      _tmp246 * _tmp281 +
      _tmp266 * (-_tmp261 * _tmp265 -
                 _tmp274 * (-_tmp246 * _tmp278 + _tmp272 * (-_tmp246 * _tmp276 - _tmp261 * _tmp277 -
                                                            _tmp268 * _tmp275))) -
      _tmp268 * _tmp280;
  _res(2, 0) =
      _tmp310 *
          (-_tmp298 * _tmp309 -
           _tmp316 * (-_tmp307 * _tmp314 + _tmp312 * (-_tmp296 * _tmp311 - _tmp309 * _tmp310))) +
      _tmp311 * _tmp317;
  _res(3, 0) =
      _tmp322 * _tmp342 +
      _tmp330 *
          (-_tmp334 * _tmp341 -
           _tmp339 * (-_tmp322 * _tmp337 + _tmp335 * (-_tmp322 * _tmp329 - _tmp330 * _tmp334)));
  _res(0, 1) = 0;
  _res(1, 1) =
      _tmp266 * (-_tmp265 * _tmp386 -
                 _tmp274 * (_tmp272 * (-_tmp267 * _tmp275 - _tmp276 * _tmp377 - _tmp277 * _tmp386) -
                            _tmp278 * _tmp377)) -
      _tmp267 * _tmp280 + _tmp281 * _tmp377;
  _res(2, 1) =
      _tmp310 *
          (-_tmp298 * _tmp391 -
           _tmp316 * (_tmp312 * (-_tmp296 * _tmp390 - _tmp310 * _tmp391) - _tmp314 * _tmp389)) +
      _tmp317 * _tmp390;
  _res(3, 1) =
      _tmp330 *
          (-_tmp339 * (_tmp335 * (-_tmp329 * _tmp392 - _tmp330 * _tmp393) - _tmp337 * _tmp392) -
           _tmp341 * _tmp393) +
      _tmp342 * _tmp392;
  _res(0, 2) = 0;
  _res(1, 2) = -_tmp279 * _tmp280;
  _res(2, 2) = 0;
  _res(3, 2) = 0;

  return _res;
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym

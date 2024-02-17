// -----------------------------------------------------------------------------
// This file was autogenerated by symforce from template:
//     function/FUNCTION.h.jinja
// Do NOT modify by hand.
// -----------------------------------------------------------------------------

#pragma once

#include <Eigen/Dense>

#include <sym/pose3.h>
#include <sym/rot3.h>

namespace sym {

/**
 * This function was autogenerated from a symbolic function. Do not modify by hand.
 *
 * Symbolic function: FK_residual_func_cost3_wrt_pb
 *
 * Args:
 *     fh1: Scalar
 *     fv1: Scalar
 *     DeltaRot: Rot3
 *     TransformationMatrix: Pose3
 *     encoder: Matrix41
 *     offset: Matrix41
 *     p_a: Matrix31
 *     p_b: Matrix31
 *     p_c: Matrix31
 *     p_d: Matrix31
 *     epsilon: Scalar
 *
 * Outputs:
 *     res: Matrix43
 */
template <typename Scalar>
Eigen::Matrix<Scalar, 4, 3> FkResidualFuncCost3WrtPb(
    const Scalar fh1, const Scalar fv1, const sym::Rot3<Scalar>& DeltaRot,
    const sym::Pose3<Scalar>& TransformationMatrix, const Eigen::Matrix<Scalar, 4, 1>& encoder,
    const Eigen::Matrix<Scalar, 4, 1>& offset, const Eigen::Matrix<Scalar, 3, 1>& p_a,
    const Eigen::Matrix<Scalar, 3, 1>& p_b, const Eigen::Matrix<Scalar, 3, 1>& p_c,
    const Eigen::Matrix<Scalar, 3, 1>& p_d, const Scalar epsilon) {
  // Total ops: 1192

  // Unused inputs
  (void)encoder;
  (void)offset;
  (void)epsilon;

  // Input arrays
  const Eigen::Matrix<Scalar, 4, 1>& _DeltaRot = DeltaRot.Data();
  const Eigen::Matrix<Scalar, 7, 1>& _TransformationMatrix = TransformationMatrix.Data();

  // Intermediate terms (387)
  const Scalar _tmp0 = Scalar(43.164000000000001) - fv1;
  const Scalar _tmp1 =
      _DeltaRot[0] * _TransformationMatrix[2] + _DeltaRot[1] * _TransformationMatrix[3] -
      _DeltaRot[2] * _TransformationMatrix[0] + _DeltaRot[3] * _TransformationMatrix[1];
  const Scalar _tmp2 =
      _DeltaRot[0] * _TransformationMatrix[3] - _DeltaRot[1] * _TransformationMatrix[2] +
      _DeltaRot[2] * _TransformationMatrix[1] + _DeltaRot[3] * _TransformationMatrix[0];
  const Scalar _tmp3 = 2 * _tmp2;
  const Scalar _tmp4 = _tmp1 * _tmp3;
  const Scalar _tmp5 =
      -_DeltaRot[0] * _TransformationMatrix[1] + _DeltaRot[1] * _TransformationMatrix[0] +
      _DeltaRot[2] * _TransformationMatrix[3] + _DeltaRot[3] * _TransformationMatrix[2];
  const Scalar _tmp6 =
      -2 * _DeltaRot[0] * _TransformationMatrix[0] - 2 * _DeltaRot[1] * _TransformationMatrix[1] -
      2 * _DeltaRot[2] * _TransformationMatrix[2] + 2 * _DeltaRot[3] * _TransformationMatrix[3];
  const Scalar _tmp7 = _tmp5 * _tmp6;
  const Scalar _tmp8 = Scalar(0.20999999999999999) * _tmp4 + Scalar(0.20999999999999999) * _tmp7;
  const Scalar _tmp9 = 2 * _tmp1 * _tmp5;
  const Scalar _tmp10 = _tmp2 * _tmp6;
  const Scalar _tmp11 = -_tmp10 + _tmp9;
  const Scalar _tmp12 = -Scalar(0.010999999999999999) * _tmp11;
  const Scalar _tmp13 = -2 * std::pow(_tmp2, Scalar(2));
  const Scalar _tmp14 = -2 * std::pow(_tmp5, Scalar(2));
  const Scalar _tmp15 = Scalar(0.20999999999999999) * _tmp13 +
                        Scalar(0.20999999999999999) * _tmp14 + Scalar(0.20999999999999999);
  const Scalar _tmp16 = _tmp12 - _tmp15;
  const Scalar _tmp17 = _tmp16 + _tmp8;
  const Scalar _tmp18 = _TransformationMatrix[5] + _tmp17;
  const Scalar _tmp19 = _tmp18 - p_b(1, 0);
  const Scalar _tmp20 = 1 - 2 * std::pow(_tmp1, Scalar(2));
  const Scalar _tmp21 = Scalar(0.20999999999999999) * _tmp14 + Scalar(0.20999999999999999) * _tmp20;
  const Scalar _tmp22 = _tmp3 * _tmp5;
  const Scalar _tmp23 = _tmp1 * _tmp6;
  const Scalar _tmp24 = _tmp22 + _tmp23;
  const Scalar _tmp25 = -Scalar(0.010999999999999999) * _tmp24;
  const Scalar _tmp26 = Scalar(0.20999999999999999) * _tmp4 - Scalar(0.20999999999999999) * _tmp7;
  const Scalar _tmp27 = _tmp25 - _tmp26;
  const Scalar _tmp28 = _tmp21 + _tmp27;
  const Scalar _tmp29 = _TransformationMatrix[4] + _tmp28;
  const Scalar _tmp30 = _tmp29 - p_b(0, 0);
  const Scalar _tmp31 = Scalar(1.0) / (_tmp30);
  const Scalar _tmp32 = _tmp19 * _tmp31;
  const Scalar _tmp33 = _tmp25 + _tmp26;
  const Scalar _tmp34 = _tmp21 + _tmp33;
  const Scalar _tmp35 = _TransformationMatrix[4] + _tmp34;
  const Scalar _tmp36 = _tmp35 - p_c(0, 0);
  const Scalar _tmp37 = std::pow(_tmp36, Scalar(2));
  const Scalar _tmp38 = _tmp12 + _tmp15;
  const Scalar _tmp39 = _tmp38 + _tmp8;
  const Scalar _tmp40 = _TransformationMatrix[5] + _tmp39;
  const Scalar _tmp41 = _tmp40 - p_c(1, 0);
  const Scalar _tmp42 = _tmp37 + std::pow(_tmp41, Scalar(2));
  const Scalar _tmp43 = std::pow(_tmp42, Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp44 = _tmp36 * _tmp43;
  const Scalar _tmp45 = _tmp41 * _tmp43;
  const Scalar _tmp46 = _tmp32 * _tmp44 - _tmp45;
  const Scalar _tmp47 = Scalar(1.0) / (_tmp46);
  const Scalar _tmp48 = _tmp44 * _tmp47;
  const Scalar _tmp49 = Scalar(1.0) * _tmp28;
  const Scalar _tmp50 = Scalar(1.0) * _tmp17;
  const Scalar _tmp51 = -_tmp50;
  const Scalar _tmp52 = Scalar(1.0) / (_tmp39 + _tmp51);
  const Scalar _tmp53 = _tmp52 * (-_tmp34 + _tmp49);
  const Scalar _tmp54 = _tmp49 + _tmp50 * _tmp53;
  const Scalar _tmp55 = 0;
  const Scalar _tmp56 = Scalar(0.20999999999999999) * _tmp22 - Scalar(0.20999999999999999) * _tmp23;
  const Scalar _tmp57 = -_tmp56;
  const Scalar _tmp58 =
      -Scalar(0.010999999999999999) * _tmp13 - Scalar(0.010999999999999999) * _tmp20;
  const Scalar _tmp59 = Scalar(0.20999999999999999) * _tmp10 + Scalar(0.20999999999999999) * _tmp9;
  const Scalar _tmp60 = _tmp58 + _tmp59;
  const Scalar _tmp61 = _tmp57 + _tmp60;
  const Scalar _tmp62 = -_tmp8;
  const Scalar _tmp63 = _tmp38 + _tmp62;
  const Scalar _tmp64 = _TransformationMatrix[5] + _tmp63;
  const Scalar _tmp65 = _tmp64 - p_d(1, 0);
  const Scalar _tmp66 = -_tmp21;
  const Scalar _tmp67 = _tmp33 + _tmp66;
  const Scalar _tmp68 = _TransformationMatrix[4] + _tmp67;
  const Scalar _tmp69 = _tmp68 - p_d(0, 0);
  const Scalar _tmp70 = std::pow(Scalar(std::pow(_tmp65, Scalar(2)) + std::pow(_tmp69, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp71 = _tmp69 * _tmp70;
  const Scalar _tmp72 = _tmp58 - _tmp59;
  const Scalar _tmp73 = _tmp56 + _tmp72;
  const Scalar _tmp74 = _tmp56 + _tmp60;
  const Scalar _tmp75 = _tmp44 * _tmp73;
  const Scalar _tmp76 = -_tmp44 * _tmp74 + _tmp75;
  const Scalar _tmp77 = _tmp65 * _tmp70;
  const Scalar _tmp78 = _tmp32 * _tmp71 - _tmp77;
  const Scalar _tmp79 = _tmp47 * _tmp78;
  const Scalar _tmp80 = _tmp31 * _tmp73;
  const Scalar _tmp81 = _tmp19 * _tmp80;
  const Scalar _tmp82 = -_tmp44 * _tmp81 + _tmp45 * _tmp74;
  const Scalar _tmp83 = _tmp61 * _tmp77 - _tmp71 * _tmp81 - _tmp79 * _tmp82;
  const Scalar _tmp84 = -_tmp53 * _tmp83 - _tmp61 * _tmp71 + _tmp71 * _tmp73 - _tmp76 * _tmp79;
  const Scalar _tmp85 = Scalar(1.0) / (_tmp84);
  const Scalar _tmp86 = _tmp78 * _tmp85;
  const Scalar _tmp87 = _tmp55 * _tmp86;
  const Scalar _tmp88 = _tmp55 * _tmp85;
  const Scalar _tmp89 = _tmp71 * _tmp88;
  const Scalar _tmp90 = _tmp0 * (-_tmp48 * _tmp87 + _tmp89);
  const Scalar _tmp91 = std::pow(_tmp30, Scalar(2));
  const Scalar _tmp92 = Scalar(1.0) / (_tmp91);
  const Scalar _tmp93 = std::pow(_tmp19, Scalar(2));
  const Scalar _tmp94 = _tmp91 + _tmp93;
  const Scalar _tmp95 = std::sqrt(_tmp94);
  const Scalar _tmp96 = _tmp92 * _tmp95;
  const Scalar _tmp97 = Scalar(1.0) / (_tmp95);
  const Scalar _tmp98 = _tmp28 * _tmp97;
  const Scalar _tmp99 = _tmp17 * _tmp97;
  const Scalar _tmp100 = _tmp19 * _tmp98 - _tmp30 * _tmp99;
  const Scalar _tmp101 = _tmp100 * _tmp95;
  const Scalar _tmp102 = _tmp101 * _tmp31;
  const Scalar _tmp103 = _tmp102 * _tmp44 - _tmp34 * _tmp45 + _tmp39 * _tmp44;
  const Scalar _tmp104 = _tmp102 * _tmp71 - _tmp103 * _tmp79 + _tmp63 * _tmp71 - _tmp67 * _tmp77;
  const Scalar _tmp105 = Scalar(1.0) / (_tmp104);
  const Scalar _tmp106 = Scalar(1.0) * _tmp105;
  const Scalar _tmp107 = _tmp106 * _tmp71;
  const Scalar _tmp108 = _tmp44 * _tmp79;
  const Scalar _tmp109 = _tmp27 + _tmp66;
  const Scalar _tmp110 = _TransformationMatrix[4] + _tmp109 - p_a(0, 0);
  const Scalar _tmp111 = _tmp16 + _tmp62;
  const Scalar _tmp112 = _TransformationMatrix[5] + _tmp111 - p_a(1, 0);
  const Scalar _tmp113 =
      std::pow(Scalar(std::pow(_tmp110, Scalar(2)) + std::pow(_tmp112, Scalar(2))),
               Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp114 = _tmp112 * _tmp113;
  const Scalar _tmp115 = _tmp110 * _tmp113;
  const Scalar _tmp116 = fh1 * (_tmp109 * _tmp114 - _tmp111 * _tmp115);
  const Scalar _tmp117 = _tmp116 * (-_tmp106 * _tmp108 + _tmp107);
  const Scalar _tmp118 = _tmp117 * _tmp97;
  const Scalar _tmp119 = _tmp19 * _tmp92;
  const Scalar _tmp120 = _tmp119 * _tmp48;
  const Scalar _tmp121 = _tmp37 / _tmp42;
  const Scalar _tmp122 = std::pow(_tmp46, Scalar(-2));
  const Scalar _tmp123 = _tmp119 * _tmp122;
  const Scalar _tmp124 = _tmp121 * _tmp123;
  const Scalar _tmp125 = _tmp106 * _tmp78;
  const Scalar _tmp126 = _tmp119 * _tmp47;
  const Scalar _tmp127 = _tmp103 * _tmp126;
  const Scalar _tmp128 = std::pow(_tmp94, Scalar(Scalar(-3) / Scalar(2)));
  const Scalar _tmp129 = _tmp128 * _tmp17;
  const Scalar _tmp130 = _tmp128 * _tmp28;
  const Scalar _tmp131 = _tmp19 * _tmp30;
  const Scalar _tmp132 = _tmp31 * _tmp95;
  const Scalar _tmp133 = _tmp132 * (-_tmp129 * _tmp91 + _tmp130 * _tmp131 + _tmp99);
  const Scalar _tmp134 = _tmp100 * _tmp97;
  const Scalar _tmp135 = _tmp101 * _tmp92;
  const Scalar _tmp136 = _tmp133 * _tmp44 - _tmp134 * _tmp44 + _tmp135 * _tmp44;
  const Scalar _tmp137 = _tmp122 * _tmp44;
  const Scalar _tmp138 = _tmp103 * _tmp137;
  const Scalar _tmp139 = _tmp119 * _tmp138;
  const Scalar _tmp140 = -_tmp127 * _tmp71 + _tmp133 * _tmp71 - _tmp134 * _tmp71 +
                         _tmp135 * _tmp71 - _tmp136 * _tmp79 + _tmp139 * _tmp78;
  const Scalar _tmp141 = std::pow(_tmp104, Scalar(-2));
  const Scalar _tmp142 = _tmp140 * _tmp141;
  const Scalar _tmp143 = Scalar(1.0) * _tmp108;
  const Scalar _tmp144 = Scalar(1.0) * _tmp71;
  const Scalar _tmp145 = _tmp116 * _tmp132;
  const Scalar _tmp146 = Scalar(1.0) * _tmp47;
  const Scalar _tmp147 = _tmp146 * _tmp82;
  const Scalar _tmp148 = -_tmp146 * _tmp76 + _tmp147 * _tmp53;
  const Scalar _tmp149 = _tmp148 * _tmp85;
  const Scalar _tmp150 = -_tmp103 * _tmp146 - _tmp104 * _tmp149;
  const Scalar _tmp151 = _tmp105 * _tmp150;
  const Scalar _tmp152 = _tmp151 * _tmp84;
  const Scalar _tmp153 = _tmp148 + _tmp152;
  const Scalar _tmp154 = -_tmp153 * _tmp86 + Scalar(1.0);
  const Scalar _tmp155 = _tmp153 * _tmp85;
  const Scalar _tmp156 = _tmp155 * _tmp71;
  const Scalar _tmp157 = _tmp114 * fh1;
  const Scalar _tmp158 = _tmp157 * (_tmp154 * _tmp48 + _tmp156);
  const Scalar _tmp159 = _tmp126 * _tmp76;
  const Scalar _tmp160 = _tmp119 * _tmp75;
  const Scalar _tmp161 = _tmp137 * _tmp82;
  const Scalar _tmp162 = _tmp119 * _tmp161;
  const Scalar _tmp163 = _tmp119 * _tmp73;
  const Scalar _tmp164 = _tmp126 * _tmp82;
  const Scalar _tmp165 = _tmp160 * _tmp79 + _tmp162 * _tmp78 - _tmp163 * _tmp71 - _tmp164 * _tmp71;
  const Scalar _tmp166 = _tmp137 * _tmp76;
  const Scalar _tmp167 = _tmp119 * _tmp166;
  const Scalar _tmp168 = -_tmp159 * _tmp71 - _tmp165 * _tmp53 + _tmp167 * _tmp78;
  const Scalar _tmp169 = _tmp32 * _tmp47;
  const Scalar _tmp170 = _tmp169 * _tmp82 + _tmp81;
  const Scalar _tmp171 = _tmp169 * _tmp76 - _tmp170 * _tmp53 - _tmp73;
  const Scalar _tmp172 = _tmp171 * _tmp85;
  const Scalar _tmp173 = -_tmp102 + _tmp103 * _tmp169 - _tmp104 * _tmp172;
  const Scalar _tmp174 = _tmp105 * _tmp173;
  const Scalar _tmp175 = _tmp168 * _tmp174;
  const Scalar _tmp176 = _tmp142 * _tmp84;
  const Scalar _tmp177 = _tmp173 * _tmp176;
  const Scalar _tmp178 = _tmp93 / [&]() {
    const Scalar base = _tmp30;
    return base * base * base;
  }();
  const Scalar _tmp179 = -_tmp161 * _tmp178 + _tmp163 + _tmp164 - _tmp178 * _tmp47 * _tmp75;
  const Scalar _tmp180 = _tmp159 - _tmp166 * _tmp178 - _tmp179 * _tmp53;
  const Scalar _tmp181 = _tmp104 * _tmp85;
  const Scalar _tmp182 = std::pow(_tmp84, Scalar(-2));
  const Scalar _tmp183 = _tmp168 * _tmp182;
  const Scalar _tmp184 = _tmp104 * _tmp171;
  const Scalar _tmp185 = _tmp105 * _tmp84;
  const Scalar _tmp186 =
      _tmp185 * (_tmp127 - _tmp133 + _tmp134 - _tmp135 + _tmp136 * _tmp169 - _tmp138 * _tmp178 -
                 _tmp140 * _tmp172 - _tmp180 * _tmp181 + _tmp183 * _tmp184);
  const Scalar _tmp187 = _tmp175 - _tmp177 + _tmp180 + _tmp186;
  const Scalar _tmp188 = _tmp71 * _tmp85;
  const Scalar _tmp189 = _tmp174 * _tmp84;
  const Scalar _tmp190 = _tmp171 + _tmp189;
  const Scalar _tmp191 = _tmp190 * _tmp85;
  const Scalar _tmp192 = _tmp191 * _tmp71;
  const Scalar _tmp193 = _tmp183 * _tmp78;
  const Scalar _tmp194 = -_tmp119 * _tmp192 - _tmp119 - _tmp187 * _tmp86 + _tmp190 * _tmp193;
  const Scalar _tmp195 = -_tmp191 * _tmp78 - _tmp32;
  const Scalar _tmp196 = _tmp121 * _tmp195;
  const Scalar _tmp197 = _tmp183 * _tmp71;
  const Scalar _tmp198 = _tmp115 * fh1;
  const Scalar _tmp199 = _tmp132 * _tmp198;
  const Scalar _tmp200 = _tmp108 * _tmp55;
  const Scalar _tmp201 = _tmp0 * _tmp132;
  const Scalar _tmp202 = _tmp198 * (_tmp192 + _tmp195 * _tmp48 + Scalar(1.0));
  const Scalar _tmp203 = _tmp202 * _tmp97;
  const Scalar _tmp204 = _tmp104 * _tmp148;
  const Scalar _tmp205 = Scalar(1.0) * _tmp162;
  const Scalar _tmp206 = _tmp146 * _tmp160;
  const Scalar _tmp207 = Scalar(1.0) * _tmp167 - _tmp205 * _tmp53 - _tmp206 * _tmp53;
  const Scalar _tmp208 = _tmp185 * (-_tmp136 * _tmp146 + Scalar(1.0) * _tmp139 - _tmp140 * _tmp149 -
                                    _tmp181 * _tmp207 + _tmp183 * _tmp204);
  const Scalar _tmp209 = _tmp150 * _tmp176;
  const Scalar _tmp210 = _tmp151 * _tmp168;
  const Scalar _tmp211 = _tmp207 + _tmp208 - _tmp209 + _tmp210;
  const Scalar _tmp212 = -_tmp119 * _tmp156 + _tmp153 * _tmp193 - _tmp211 * _tmp86;
  const Scalar _tmp213 = _tmp132 * _tmp157;
  const Scalar _tmp214 = _tmp90 * _tmp97;
  const Scalar _tmp215 = _tmp158 * _tmp97;
  const Scalar _tmp216 =
      -_tmp117 * _tmp96 + _tmp118 -
      _tmp145 * (-_tmp107 * _tmp120 + _tmp124 * _tmp125 + _tmp142 * _tmp143 - _tmp142 * _tmp144) -
      _tmp158 * _tmp96 -
      _tmp199 * (-_tmp123 * _tmp196 + _tmp187 * _tmp188 - _tmp190 * _tmp197 + _tmp194 * _tmp48) -
      _tmp201 * (-_tmp120 * _tmp89 + _tmp124 * _tmp87 + _tmp183 * _tmp200 - _tmp197 * _tmp55) -
      _tmp202 * _tmp96 + _tmp203 -
      _tmp213 * (-_tmp124 * _tmp154 - _tmp153 * _tmp197 + _tmp188 * _tmp211 + _tmp212 * _tmp48) +
      _tmp214 + _tmp215 - _tmp90 * _tmp96;
  const Scalar _tmp217 = -_tmp29 + p_b(0, 0);
  const Scalar _tmp218 = -_tmp18 + p_b(1, 0);
  const Scalar _tmp219 =
      std::sqrt(Scalar(std::pow(_tmp217, Scalar(2)) + std::pow(_tmp218, Scalar(2))));
  const Scalar _tmp220 =
      -_tmp117 * _tmp132 - _tmp132 * _tmp158 - _tmp132 * _tmp202 - _tmp132 * _tmp90;
  const Scalar _tmp221 = Scalar(1.0) / (_tmp220);
  const Scalar _tmp222 = _tmp51 + _tmp63;
  const Scalar _tmp223 = _tmp222 * _tmp53;
  const Scalar _tmp224 = Scalar(1.0) / (-_tmp223 + _tmp49 - _tmp67);
  const Scalar _tmp225 = Scalar(1.0) * _tmp224;
  const Scalar _tmp226 = _tmp185 * _tmp225;
  const Scalar _tmp227 = -_tmp106 * _tmp83 + _tmp222 * _tmp226;
  const Scalar _tmp228 = Scalar(1.0) * _tmp52;
  const Scalar _tmp229 = Scalar(1.0) * _tmp116;
  const Scalar _tmp230 = _tmp83 * _tmp85;
  const Scalar _tmp231 = _tmp222 * _tmp224;
  const Scalar _tmp232 = _tmp170 + _tmp189 * _tmp231 - _tmp190 * _tmp230;
  const Scalar _tmp233 = Scalar(1.0) * _tmp198;
  const Scalar _tmp234 = fh1 * (_tmp57 + _tmp72);
  const Scalar _tmp235 = -Scalar(5.1796800000000003) * _tmp11 - _tmp111 * fv1 - _tmp114 * _tmp234;
  const Scalar _tmp236 = _tmp225 * _tmp53;
  const Scalar _tmp237 = _tmp223 * _tmp225 + Scalar(1.0);
  const Scalar _tmp238 = -_tmp147 + _tmp152 * _tmp231 - _tmp153 * _tmp230;
  const Scalar _tmp239 = Scalar(1.0) * _tmp157;
  const Scalar _tmp240 = _tmp224 * _tmp54;
  const Scalar _tmp241 = -_tmp222 * _tmp240 - _tmp230 * _tmp55 + _tmp51;
  const Scalar _tmp242 = _tmp222 * _tmp52;
  const Scalar _tmp243 = _tmp109 * fv1 + _tmp115 * _tmp234 + Scalar(5.1796800000000003) * _tmp24;
  const Scalar _tmp244 =
      Scalar(1.0) * _tmp0 * (-_tmp225 * _tmp54 - _tmp228 * _tmp241 + Scalar(1.0)) +
      _tmp229 * (_tmp226 - _tmp227 * _tmp228) + _tmp233 * (_tmp189 * _tmp225 - _tmp228 * _tmp232) +
      Scalar(1.0) * _tmp235 * (-_tmp228 * _tmp237 + _tmp236) +
      _tmp239 * (_tmp152 * _tmp225 - _tmp228 * _tmp238) +
      Scalar(1.0) * _tmp243 * (_tmp225 * _tmp242 - _tmp225);
  const Scalar _tmp245 = std::asinh(_tmp221 * _tmp244);
  const Scalar _tmp246 = Scalar(9.6622558468725703) * _tmp245;
  const Scalar _tmp247 = -_tmp219 - _tmp220 * _tmp246;
  const Scalar _tmp248 = Scalar(0.1034955) * _tmp221;
  const Scalar _tmp249 = _tmp247 * _tmp248;
  const Scalar _tmp250 = Scalar(1.0) * _tmp245;
  const Scalar _tmp251 = -Scalar(9.6622558468725703) * std::sinh(_tmp249) -
                         Scalar(9.6622558468725703) * std::sinh(_tmp250);
  const Scalar _tmp252 = std::pow(_tmp220, Scalar(-2));
  const Scalar _tmp253 = Scalar(0.1034955) * _tmp247 * _tmp252;
  const Scalar _tmp254 = Scalar(1.0) / (_tmp219);
  const Scalar _tmp255 = _tmp244 * _tmp252;
  const Scalar _tmp256 = _tmp165 * _tmp85;
  const Scalar _tmp257 = _tmp55 * _tmp83;
  const Scalar _tmp258 = _tmp0 * _tmp52;
  const Scalar _tmp259 = _tmp258 * (_tmp183 * _tmp257 - _tmp256 * _tmp55);
  const Scalar _tmp260 = _tmp105 * _tmp225;
  const Scalar _tmp261 = _tmp168 * _tmp260;
  const Scalar _tmp262 = _tmp176 * _tmp225;
  const Scalar _tmp263 = Scalar(1.0) * _tmp83;
  const Scalar _tmp264 =
      -_tmp106 * _tmp165 + _tmp142 * _tmp263 + _tmp222 * _tmp261 - _tmp222 * _tmp262;
  const Scalar _tmp265 = _tmp153 * _tmp83;
  const Scalar _tmp266 = -_tmp153 * _tmp256 + _tmp183 * _tmp265 + _tmp205 + _tmp206 +
                         _tmp208 * _tmp231 - _tmp209 * _tmp231 + _tmp210 * _tmp231 -
                         _tmp211 * _tmp230;
  const Scalar _tmp267 = _tmp190 * _tmp83;
  const Scalar _tmp268 = _tmp175 * _tmp231 - _tmp177 * _tmp231 + _tmp179 + _tmp183 * _tmp267 +
                         _tmp186 * _tmp231 - _tmp187 * _tmp230 - _tmp190 * _tmp256;
  const Scalar _tmp269 =
      -_tmp216 * _tmp255 + _tmp221 * (_tmp229 * (-_tmp228 * _tmp264 + _tmp261 - _tmp262) +
                                      _tmp233 * (_tmp175 * _tmp225 - _tmp177 * _tmp225 +
                                                 _tmp186 * _tmp225 - _tmp228 * _tmp268) +
                                      _tmp239 * (_tmp208 * _tmp225 - _tmp209 * _tmp225 +
                                                 _tmp210 * _tmp225 - _tmp228 * _tmp266) -
                                      Scalar(1.0) * _tmp259);
  const Scalar _tmp270 =
      std::pow(Scalar(std::pow(_tmp244, Scalar(2)) * _tmp252 + 1), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp271 = Scalar(9.6622558468725703) * _tmp220;
  const Scalar _tmp272 = _tmp270 * _tmp271;
  const Scalar _tmp273 = std::cosh(_tmp249);
  const Scalar _tmp274 = Scalar(1.0) * _tmp270 * std::cosh(_tmp250);
  const Scalar _tmp275 = _tmp0 * _tmp55;
  const Scalar _tmp276 = _tmp275 * _tmp86;
  const Scalar _tmp277 = _tmp106 * _tmp116;
  const Scalar _tmp278 = _tmp198 * _tmp47;
  const Scalar _tmp279 = _tmp157 * _tmp47;
  const Scalar _tmp280 =
      _tmp154 * _tmp279 + _tmp195 * _tmp278 - _tmp276 * _tmp47 - _tmp277 * _tmp79;
  const Scalar _tmp281 = Scalar(1.0) / (_tmp280);
  const Scalar _tmp282 = _tmp116 * _tmp52;
  const Scalar _tmp283 = _tmp225 * _tmp243;
  const Scalar _tmp284 = _tmp198 * _tmp52;
  const Scalar _tmp285 = _tmp157 * _tmp52;
  const Scalar _tmp286 = _tmp227 * _tmp282 + _tmp232 * _tmp284 + _tmp235 * _tmp237 * _tmp52 +
                         _tmp238 * _tmp285 + _tmp241 * _tmp258 - _tmp242 * _tmp283;
  const Scalar _tmp287 = std::asinh(_tmp281 * _tmp286);
  const Scalar _tmp288 = Scalar(1.0) * _tmp287;
  const Scalar _tmp289 = Scalar(1.0) * std::cosh(_tmp288);
  const Scalar _tmp290 = std::pow(_tmp280, Scalar(-2));
  const Scalar _tmp291 =
      std::pow(Scalar(std::pow(_tmp286, Scalar(2)) * _tmp290 + 1), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp292 = _tmp119 * _tmp137;
  const Scalar _tmp293 = _tmp0 * _tmp88;
  const Scalar _tmp294 = _tmp126 * _tmp71;
  const Scalar _tmp295 = _tmp183 * _tmp275;
  const Scalar _tmp296 = _tmp142 * _tmp229;
  const Scalar _tmp297 = _tmp195 * _tmp198;
  const Scalar _tmp298 = -_tmp154 * _tmp157 * _tmp292 + _tmp194 * _tmp278 + _tmp212 * _tmp279 +
                         _tmp276 * _tmp292 + _tmp277 * _tmp292 * _tmp78 - _tmp277 * _tmp294 -
                         _tmp292 * _tmp297 - _tmp293 * _tmp294 + _tmp295 * _tmp79 +
                         _tmp296 * _tmp79;
  const Scalar _tmp299 = _tmp286 * _tmp290;
  const Scalar _tmp300 =
      _tmp291 * (_tmp281 * (_tmp259 + _tmp264 * _tmp282 + _tmp266 * _tmp285 + _tmp268 * _tmp284) -
                 _tmp298 * _tmp299);
  const Scalar _tmp301 = Scalar(9.6622558468725703) * _tmp287;
  const Scalar _tmp302 =
      -_tmp280 * _tmp301 - std::sqrt(Scalar(std::pow(Scalar(-_tmp35 + p_c(0, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp40 + p_c(1, 0)), Scalar(2))));
  const Scalar _tmp303 = Scalar(0.1034955) * _tmp290 * _tmp302;
  const Scalar _tmp304 = Scalar(9.6622558468725703) * _tmp280;
  const Scalar _tmp305 = Scalar(0.1034955) * _tmp281;
  const Scalar _tmp306 = _tmp302 * _tmp305;
  const Scalar _tmp307 = std::cosh(_tmp306);
  const Scalar _tmp308 = -Scalar(9.6622558468725703) * std::sinh(_tmp288) -
                         Scalar(9.6622558468725703) * std::sinh(_tmp306);
  const Scalar _tmp309 = _tmp155 * _tmp157 + _tmp191 * _tmp198 + _tmp277 + _tmp293;
  const Scalar _tmp310 = Scalar(1.0) / (_tmp309);
  const Scalar _tmp311 = _tmp157 * _tmp224;
  const Scalar _tmp312 = _tmp198 * _tmp224;
  const Scalar _tmp313 = _tmp0 * _tmp240 - _tmp116 * _tmp226 - _tmp152 * _tmp311 -
                         _tmp189 * _tmp312 - _tmp235 * _tmp236 + _tmp283;
  const Scalar _tmp314 = std::asinh(_tmp310 * _tmp313);
  const Scalar _tmp315 = Scalar(9.6622558468725703) * _tmp309;
  const Scalar _tmp316 =
      -_tmp314 * _tmp315 - std::sqrt(Scalar(std::pow(Scalar(-_tmp64 + p_d(1, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp68 + p_d(0, 0)), Scalar(2))));
  const Scalar _tmp317 = Scalar(0.1034955) * _tmp310;
  const Scalar _tmp318 = _tmp316 * _tmp317;
  const Scalar _tmp319 = Scalar(1.0) * _tmp314;
  const Scalar _tmp320 = -std::sinh(_tmp318) - std::sinh(_tmp319);
  const Scalar _tmp321 = _tmp153 * _tmp157;
  const Scalar _tmp322 = _tmp190 * _tmp198;
  const Scalar _tmp323 = _tmp198 * _tmp85;
  const Scalar _tmp324 = _tmp157 * _tmp85;
  const Scalar _tmp325 = -_tmp183 * _tmp321 - _tmp183 * _tmp322 + _tmp187 * _tmp323 +
                         _tmp211 * _tmp324 - _tmp295 - _tmp296;
  const Scalar _tmp326 = Scalar(9.6622558468725703) * _tmp325;
  const Scalar _tmp327 = std::pow(_tmp309, Scalar(-2));
  const Scalar _tmp328 = _tmp313 * _tmp327;
  const Scalar _tmp329 =
      std::pow(Scalar(std::pow(_tmp313, Scalar(2)) * _tmp327 + 1), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp330 =
      _tmp329 *
      (_tmp310 * (-_tmp116 * _tmp261 + _tmp116 * _tmp262 - _tmp175 * _tmp312 + _tmp177 * _tmp312 -
                  _tmp186 * _tmp312 - _tmp208 * _tmp311 + _tmp209 * _tmp311 - _tmp210 * _tmp311) -
       _tmp325 * _tmp328);
  const Scalar _tmp331 = Scalar(1.0) * std::cosh(_tmp319);
  const Scalar _tmp332 = Scalar(0.1034955) * _tmp316 * _tmp327;
  const Scalar _tmp333 = std::cosh(_tmp318);
  const Scalar _tmp334 = _tmp132 * (-_tmp129 * _tmp131 + _tmp130 * _tmp93 - _tmp98);
  const Scalar _tmp335 = _tmp134 * _tmp32;
  const Scalar _tmp336 = _tmp334 * _tmp44 - _tmp335 * _tmp44;
  const Scalar _tmp337 = _tmp122 * _tmp31;
  const Scalar _tmp338 = _tmp337 * _tmp44;
  const Scalar _tmp339 = Scalar(1.0) * _tmp338;
  const Scalar _tmp340 = _tmp339 * _tmp82;
  const Scalar _tmp341 = _tmp146 * _tmp44 * _tmp80;
  const Scalar _tmp342 = -_tmp339 * _tmp76 + _tmp340 * _tmp53 + _tmp341 * _tmp53;
  const Scalar _tmp343 = _tmp338 * _tmp78;
  const Scalar _tmp344 = _tmp31 * _tmp47;
  const Scalar _tmp345 = _tmp344 * _tmp76;
  const Scalar _tmp346 = _tmp344 * _tmp82;
  const Scalar _tmp347 = -_tmp108 * _tmp80 - _tmp343 * _tmp82 + _tmp346 * _tmp71 + _tmp71 * _tmp80;
  const Scalar _tmp348 = -_tmp343 * _tmp76 + _tmp345 * _tmp71 - _tmp347 * _tmp53;
  const Scalar _tmp349 = _tmp182 * _tmp348;
  const Scalar _tmp350 = _tmp103 * _tmp344;
  const Scalar _tmp351 = -_tmp103 * _tmp343 + _tmp334 * _tmp71 - _tmp335 * _tmp71 -
                         _tmp336 * _tmp79 + _tmp350 * _tmp71;
  const Scalar _tmp352 = _tmp185 * (-_tmp103 * _tmp339 - _tmp146 * _tmp336 - _tmp149 * _tmp351 -
                                    _tmp181 * _tmp342 + _tmp204 * _tmp349);
  const Scalar _tmp353 = _tmp141 * _tmp351;
  const Scalar _tmp354 = _tmp353 * _tmp84;
  const Scalar _tmp355 = _tmp150 * _tmp354;
  const Scalar _tmp356 = _tmp151 * _tmp348;
  const Scalar _tmp357 = _tmp342 + _tmp352 - _tmp355 + _tmp356;
  const Scalar _tmp358 = -_tmp155 * _tmp347 - _tmp230 * _tmp357 + _tmp231 * _tmp352 -
                         _tmp231 * _tmp355 + _tmp231 * _tmp356 + _tmp265 * _tmp349 - _tmp340 -
                         _tmp341;
  const Scalar _tmp359 = _tmp126 * _tmp75 + _tmp162 - _tmp346 - _tmp80;
  const Scalar _tmp360 = _tmp167 - _tmp345 - _tmp359 * _tmp53;
  const Scalar _tmp361 =
      _tmp185 * (_tmp139 + _tmp169 * _tmp336 - _tmp172 * _tmp351 - _tmp181 * _tmp360 +
                 _tmp184 * _tmp349 - _tmp334 + _tmp335 - _tmp350);
  const Scalar _tmp362 = _tmp173 * _tmp354;
  const Scalar _tmp363 = _tmp174 * _tmp348;
  const Scalar _tmp364 = _tmp360 + _tmp361 - _tmp362 + _tmp363;
  const Scalar _tmp365 = -_tmp191 * _tmp347 - _tmp230 * _tmp364 + _tmp231 * _tmp361 -
                         _tmp231 * _tmp362 + _tmp231 * _tmp363 + _tmp267 * _tmp349 + _tmp359;
  const Scalar _tmp366 = _tmp258 * (_tmp257 * _tmp349 - _tmp347 * _tmp88);
  const Scalar _tmp367 = _tmp260 * _tmp348;
  const Scalar _tmp368 = _tmp225 * _tmp354;
  const Scalar _tmp369 =
      -_tmp106 * _tmp347 + _tmp222 * _tmp367 - _tmp222 * _tmp368 + _tmp263 * _tmp353;
  const Scalar _tmp370 = _tmp349 * _tmp71;
  const Scalar _tmp371 = _tmp349 * _tmp78;
  const Scalar _tmp372 = _tmp190 * _tmp371 + _tmp192 * _tmp31 + _tmp31 - _tmp364 * _tmp86;
  const Scalar _tmp373 = _tmp121 * _tmp337;
  const Scalar _tmp374 = _tmp344 * _tmp44;
  const Scalar _tmp375 = _tmp153 * _tmp371 + _tmp156 * _tmp31 - _tmp357 * _tmp86;
  const Scalar _tmp376 = _tmp154 * _tmp337;
  const Scalar _tmp377 =
      _tmp118 * _tmp32 -
      _tmp145 * (_tmp107 * _tmp374 - _tmp125 * _tmp373 + _tmp143 * _tmp353 - _tmp144 * _tmp353) -
      _tmp199 * (_tmp188 * _tmp364 - _tmp190 * _tmp370 + _tmp196 * _tmp337 + _tmp372 * _tmp48) -
      _tmp201 * (_tmp200 * _tmp349 - _tmp370 * _tmp55 - _tmp373 * _tmp87 + _tmp374 * _tmp89) +
      _tmp203 * _tmp32 -
      _tmp213 * (_tmp121 * _tmp376 - _tmp153 * _tmp370 + _tmp188 * _tmp357 + _tmp375 * _tmp48) +
      _tmp214 * _tmp32 + _tmp215 * _tmp32;
  const Scalar _tmp378 = _tmp221 * (_tmp229 * (-_tmp228 * _tmp369 + _tmp367 - _tmp368) +
                                    _tmp233 * (_tmp225 * _tmp361 - _tmp225 * _tmp362 +
                                               _tmp225 * _tmp363 - _tmp228 * _tmp365) +
                                    _tmp239 * (_tmp225 * _tmp352 - _tmp225 * _tmp355 +
                                               _tmp225 * _tmp356 - _tmp228 * _tmp358) -
                                    Scalar(1.0) * _tmp366) -
                         _tmp255 * _tmp377;
  const Scalar _tmp379 = _tmp275 * _tmp349;
  const Scalar _tmp380 = _tmp344 * _tmp71;
  const Scalar _tmp381 = _tmp229 * _tmp353;
  const Scalar _tmp382 = _tmp157 * _tmp376 * _tmp44 - _tmp276 * _tmp338 - _tmp277 * _tmp343 +
                         _tmp277 * _tmp380 + _tmp278 * _tmp372 + _tmp279 * _tmp375 +
                         _tmp293 * _tmp380 + _tmp297 * _tmp338 + _tmp379 * _tmp79 +
                         _tmp381 * _tmp79;
  const Scalar _tmp383 =
      _tmp291 * (_tmp281 * (_tmp282 * _tmp369 + _tmp284 * _tmp365 + _tmp285 * _tmp358 + _tmp366) -
                 _tmp299 * _tmp382);
  const Scalar _tmp384 = -_tmp321 * _tmp349 - _tmp322 * _tmp349 + _tmp323 * _tmp364 +
                         _tmp324 * _tmp357 - _tmp379 - _tmp381;
  const Scalar _tmp385 = Scalar(9.6622558468725703) * _tmp384;
  const Scalar _tmp386 =
      _tmp329 *
      (_tmp310 * (-_tmp116 * _tmp367 + _tmp116 * _tmp368 - _tmp311 * _tmp352 + _tmp311 * _tmp355 -
                  _tmp311 * _tmp356 - _tmp312 * _tmp361 + _tmp312 * _tmp362 - _tmp312 * _tmp363) -
       _tmp328 * _tmp384);

  // Output terms (1)
  Eigen::Matrix<Scalar, 4, 3> _res;

  _res(0, 0) = 0;
  _res(1, 0) =
      _tmp216 * _tmp251 +
      _tmp271 * (-_tmp269 * _tmp274 -
                 _tmp273 * (-_tmp216 * _tmp253 + _tmp248 * (-_tmp216 * _tmp246 - _tmp217 * _tmp254 -
                                                            _tmp269 * _tmp272)));
  _res(2, 0) =
      _tmp298 * _tmp308 +
      _tmp304 *
          (-_tmp289 * _tmp300 -
           _tmp307 * (-_tmp298 * _tmp303 + _tmp305 * (-_tmp298 * _tmp301 - _tmp300 * _tmp304)));
  _res(3, 0) =
      _tmp315 *
          (-_tmp330 * _tmp331 -
           _tmp333 * (_tmp317 * (-_tmp314 * _tmp326 - _tmp315 * _tmp330) - _tmp325 * _tmp332)) +
      _tmp320 * _tmp326;
  _res(0, 1) = 0;
  _res(1, 1) =
      _tmp251 * _tmp377 +
      _tmp271 *
          (-_tmp273 * (_tmp248 * (-_tmp218 * _tmp254 - _tmp246 * _tmp377 - _tmp272 * _tmp378) -
                       _tmp253 * _tmp377) -
           _tmp274 * _tmp378);
  _res(2, 1) =
      _tmp304 *
          (-_tmp289 * _tmp383 -
           _tmp307 * (-_tmp303 * _tmp382 + _tmp305 * (-_tmp301 * _tmp382 - _tmp304 * _tmp383))) +
      _tmp308 * _tmp382;
  _res(3, 1) =
      _tmp315 *
          (-_tmp331 * _tmp386 -
           _tmp333 * (_tmp317 * (-_tmp314 * _tmp385 - _tmp315 * _tmp386) - _tmp332 * _tmp384)) +
      _tmp320 * _tmp385;
  _res(0, 2) = 0;
  _res(1, 2) = 0;
  _res(2, 2) = 0;
  _res(3, 2) = 0;

  return _res;
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym

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
 * Symbolic function: FK_residual_func_cost1_wrt_pa
 *
 * Args:
 *     fh1: Scalar
 *     fv1: Scalar
 *     DeltaRot: Rot3
 *     TransformationMatrix: Pose3
 *     encoder: Matrix41
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
Eigen::Matrix<Scalar, 4, 3> FkResidualFuncCost1WrtPa(
    const Scalar fh1, const Scalar fv1, const sym::Rot3<Scalar>& DeltaRot,
    const sym::Pose3<Scalar>& TransformationMatrix, const Eigen::Matrix<Scalar, 4, 1>& encoder,
    const Eigen::Matrix<Scalar, 3, 1>& p_a, const Eigen::Matrix<Scalar, 3, 1>& p_b,
    const Eigen::Matrix<Scalar, 3, 1>& p_c, const Eigen::Matrix<Scalar, 3, 1>& p_d,
    const Scalar epsilon) {
  // Total ops: 777

  // Unused inputs
  (void)encoder;
  (void)epsilon;

  // Input arrays
  const Eigen::Matrix<Scalar, 4, 1>& _DeltaRot = DeltaRot.Data();
  const Eigen::Matrix<Scalar, 7, 1>& _TransformationMatrix = TransformationMatrix.Data();

  // Intermediate terms (273)
  const Scalar _tmp0 =
      _DeltaRot[0] * _TransformationMatrix[3] - _DeltaRot[1] * _TransformationMatrix[2] +
      _DeltaRot[2] * _TransformationMatrix[1] + _DeltaRot[3] * _TransformationMatrix[0];
  const Scalar _tmp1 =
      _DeltaRot[0] * _TransformationMatrix[2] + _DeltaRot[1] * _TransformationMatrix[3] -
      _DeltaRot[2] * _TransformationMatrix[0] + _DeltaRot[3] * _TransformationMatrix[1];
  const Scalar _tmp2 = 2 * _tmp0 * _tmp1;
  const Scalar _tmp3 =
      -_DeltaRot[0] * _TransformationMatrix[1] + _DeltaRot[1] * _TransformationMatrix[0] +
      _DeltaRot[2] * _TransformationMatrix[3] + _DeltaRot[3] * _TransformationMatrix[2];
  const Scalar _tmp4 =
      -2 * _DeltaRot[0] * _TransformationMatrix[0] - 2 * _DeltaRot[1] * _TransformationMatrix[1] -
      2 * _DeltaRot[2] * _TransformationMatrix[2] + 2 * _DeltaRot[3] * _TransformationMatrix[3];
  const Scalar _tmp5 = _tmp3 * _tmp4;
  const Scalar _tmp6 = Scalar(0.20999999999999999) * _tmp2 - Scalar(0.20999999999999999) * _tmp5;
  const Scalar _tmp7 = -_tmp6;
  const Scalar _tmp8 = 2 * _tmp3;
  const Scalar _tmp9 = _tmp0 * _tmp8;
  const Scalar _tmp10 = _tmp1 * _tmp4;
  const Scalar _tmp11 = _tmp10 + _tmp9;
  const Scalar _tmp12 = -Scalar(0.010999999999999999) * _tmp11;
  const Scalar _tmp13 = -2 * std::pow(_tmp3, Scalar(2));
  const Scalar _tmp14 = 1 - 2 * std::pow(_tmp1, Scalar(2));
  const Scalar _tmp15 = Scalar(0.20999999999999999) * _tmp13 + Scalar(0.20999999999999999) * _tmp14;
  const Scalar _tmp16 = _tmp12 - _tmp15;
  const Scalar _tmp17 = _tmp16 + _tmp7;
  const Scalar _tmp18 = _TransformationMatrix[4] + _tmp17;
  const Scalar _tmp19 = -_tmp18 + p_a(0, 0);
  const Scalar _tmp20 = Scalar(1.0) / (fh1);
  const Scalar _tmp21 = Scalar(0.20999999999999999) * _tmp2 + Scalar(0.20999999999999999) * _tmp5;
  const Scalar _tmp22 = -_tmp21;
  const Scalar _tmp23 = _tmp1 * _tmp8;
  const Scalar _tmp24 = _tmp0 * _tmp4;
  const Scalar _tmp25 = _tmp23 - _tmp24;
  const Scalar _tmp26 = -Scalar(0.010999999999999999) * _tmp25;
  const Scalar _tmp27 = -2 * std::pow(_tmp0, Scalar(2));
  const Scalar _tmp28 = Scalar(0.20999999999999999) * _tmp13 +
                        Scalar(0.20999999999999999) * _tmp27 + Scalar(0.20999999999999999);
  const Scalar _tmp29 = _tmp26 - _tmp28;
  const Scalar _tmp30 = _tmp22 + _tmp29;
  const Scalar _tmp31 = _TransformationMatrix[5] + _tmp30;
  const Scalar _tmp32 = -_tmp31 + p_a(1, 0);
  const Scalar _tmp33 =
      std::sqrt(Scalar(std::pow(_tmp19, Scalar(2)) + std::pow(_tmp32, Scalar(2))));
  const Scalar _tmp34 =
      Scalar(1.0) *
      std::sinh(Scalar(0.1034955) * _tmp20 *
                (-_tmp33 - Scalar(9.6622558468725703) * fh1 * std::asinh(_tmp20 * fv1))) /
      _tmp33;
  const Scalar _tmp35 = Scalar(43.164000000000001) - fv1;
  const Scalar _tmp36 = _tmp12 + _tmp15;
  const Scalar _tmp37 = _tmp36 + _tmp6;
  const Scalar _tmp38 = _TransformationMatrix[4] + _tmp37;
  const Scalar _tmp39 = _tmp38 - p_c(0, 0);
  const Scalar _tmp40 = _tmp26 + _tmp28;
  const Scalar _tmp41 = _tmp21 + _tmp40;
  const Scalar _tmp42 = _TransformationMatrix[5] + _tmp41;
  const Scalar _tmp43 = _tmp42 - p_c(1, 0);
  const Scalar _tmp44 = std::pow(Scalar(std::pow(_tmp39, Scalar(2)) + std::pow(_tmp43, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp45 = _tmp39 * _tmp44;
  const Scalar _tmp46 = _tmp36 + _tmp7;
  const Scalar _tmp47 = Scalar(1.0) * _tmp46;
  const Scalar _tmp48 = _tmp21 + _tmp29;
  const Scalar _tmp49 = Scalar(1.0) * _tmp48;
  const Scalar _tmp50 = -_tmp49;
  const Scalar _tmp51 = Scalar(1.0) / (_tmp41 + _tmp50);
  const Scalar _tmp52 = -_tmp37 + _tmp47;
  const Scalar _tmp53 = _tmp51 * _tmp52;
  const Scalar _tmp54 = _tmp47 + _tmp49 * _tmp53;
  const Scalar _tmp55 = 0;
  const Scalar _tmp56 = _TransformationMatrix[5] + _tmp48;
  const Scalar _tmp57 = _tmp56 - p_b(1, 0);
  const Scalar _tmp58 = _TransformationMatrix[4] + _tmp46;
  const Scalar _tmp59 = _tmp58 - p_b(0, 0);
  const Scalar _tmp60 = Scalar(1.0) / (_tmp59);
  const Scalar _tmp61 = _tmp57 * _tmp60;
  const Scalar _tmp62 = _tmp43 * _tmp44;
  const Scalar _tmp63 = Scalar(1.0) / (_tmp45 * _tmp61 - _tmp62);
  const Scalar _tmp64 = Scalar(0.20999999999999999) * _tmp23 + Scalar(0.20999999999999999) * _tmp24;
  const Scalar _tmp65 =
      -Scalar(0.010999999999999999) * _tmp14 - Scalar(0.010999999999999999) * _tmp27;
  const Scalar _tmp66 = -Scalar(0.20999999999999999) * _tmp10 + Scalar(0.20999999999999999) * _tmp9;
  const Scalar _tmp67 = _tmp65 - _tmp66;
  const Scalar _tmp68 = _tmp64 + _tmp67;
  const Scalar _tmp69 = _tmp22 + _tmp40;
  const Scalar _tmp70 = _TransformationMatrix[5] + _tmp69;
  const Scalar _tmp71 = _tmp70 - p_d(1, 0);
  const Scalar _tmp72 = _tmp16 + _tmp6;
  const Scalar _tmp73 = _TransformationMatrix[4] + _tmp72;
  const Scalar _tmp74 = _tmp73 - p_d(0, 0);
  const Scalar _tmp75 = std::pow(Scalar(std::pow(_tmp71, Scalar(2)) + std::pow(_tmp74, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp76 = _tmp74 * _tmp75;
  const Scalar _tmp77 = -_tmp64;
  const Scalar _tmp78 = _tmp65 + _tmp66;
  const Scalar _tmp79 = _tmp77 + _tmp78;
  const Scalar _tmp80 = _tmp64 + _tmp78;
  const Scalar _tmp81 = _tmp45 * _tmp79 - _tmp45 * _tmp80;
  const Scalar _tmp82 = _tmp71 * _tmp75;
  const Scalar _tmp83 = _tmp61 * _tmp76 - _tmp82;
  const Scalar _tmp84 = _tmp63 * _tmp83;
  const Scalar _tmp85 = _tmp61 * _tmp79;
  const Scalar _tmp86 = _tmp63 * (-_tmp45 * _tmp85 + _tmp62 * _tmp80);
  const Scalar _tmp87 = _tmp68 * _tmp82 - _tmp76 * _tmp85 - _tmp83 * _tmp86;
  const Scalar _tmp88 = -_tmp53 * _tmp87 - _tmp68 * _tmp76 + _tmp76 * _tmp79 - _tmp81 * _tmp84;
  const Scalar _tmp89 = Scalar(1.0) / (_tmp88);
  const Scalar _tmp90 = _tmp83 * _tmp89;
  const Scalar _tmp91 = _tmp55 * _tmp63 * _tmp90;
  const Scalar _tmp92 = _tmp55 * _tmp89;
  const Scalar _tmp93 =
      std::sqrt(Scalar(std::pow(_tmp57, Scalar(2)) + std::pow(_tmp59, Scalar(2))));
  const Scalar _tmp94 = _tmp60 * _tmp93;
  const Scalar _tmp95 = _tmp31 - p_a(1, 0);
  const Scalar _tmp96 = Scalar(1.0) * _tmp63;
  const Scalar _tmp97 = Scalar(1.0) * _tmp51;
  const Scalar _tmp98 = _tmp52 * _tmp86 * _tmp97 - _tmp81 * _tmp96;
  const Scalar _tmp99 = Scalar(1.0) / (_tmp93);
  const Scalar _tmp100 = _tmp94 * (_tmp46 * _tmp57 * _tmp99 - _tmp48 * _tmp59 * _tmp99);
  const Scalar _tmp101 = _tmp100 * _tmp45 - _tmp37 * _tmp62 + _tmp41 * _tmp45;
  const Scalar _tmp102 = _tmp100 * _tmp76 - _tmp101 * _tmp84 + _tmp69 * _tmp76 - _tmp72 * _tmp82;
  const Scalar _tmp103 = _tmp102 * _tmp89;
  const Scalar _tmp104 = Scalar(1.0) / (_tmp102);
  const Scalar _tmp105 = _tmp104 * _tmp88;
  const Scalar _tmp106 = _tmp105 * (-_tmp101 * _tmp96 - _tmp103 * _tmp98);
  const Scalar _tmp107 = _tmp106 + _tmp98;
  const Scalar _tmp108 = _tmp63 * (-_tmp107 * _tmp90 + Scalar(1.0));
  const Scalar _tmp109 = _tmp107 * _tmp89;
  const Scalar _tmp110 = _tmp108 * _tmp45 + _tmp109 * _tmp76;
  const Scalar _tmp111 = _tmp18 - p_a(0, 0);
  const Scalar _tmp112 = std::pow(_tmp111, Scalar(2));
  const Scalar _tmp113 = std::pow(_tmp95, Scalar(2));
  const Scalar _tmp114 = _tmp112 + _tmp113;
  const Scalar _tmp115 = std::pow(_tmp114, Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp116 = _tmp115 * fh1;
  const Scalar _tmp117 = _tmp116 * _tmp94;
  const Scalar _tmp118 = _tmp110 * _tmp117;
  const Scalar _tmp119 = _tmp61 * _tmp63;
  const Scalar _tmp120 = _tmp61 * _tmp86 + _tmp85;
  const Scalar _tmp121 = _tmp119 * _tmp81 - _tmp120 * _tmp53 - _tmp79;
  const Scalar _tmp122 = _tmp105 * (-_tmp100 + _tmp101 * _tmp119 - _tmp103 * _tmp121);
  const Scalar _tmp123 = _tmp121 + _tmp122;
  const Scalar _tmp124 = _tmp63 * (-_tmp123 * _tmp90 - _tmp61);
  const Scalar _tmp125 = _tmp123 * _tmp89;
  const Scalar _tmp126 = _tmp124 * _tmp45 + _tmp125 * _tmp76 + Scalar(1.0);
  const Scalar _tmp127 = _tmp117 * _tmp126;
  const Scalar _tmp128 = _tmp115 * _tmp17;
  const Scalar _tmp129 = _tmp115 * _tmp30;
  const Scalar _tmp130 = fh1 * (-_tmp111 * _tmp129 + _tmp128 * _tmp95);
  const Scalar _tmp131 = Scalar(1.0) * _tmp104;
  const Scalar _tmp132 = _tmp94 * (-_tmp131 * _tmp45 * _tmp84 + _tmp131 * _tmp76);
  const Scalar _tmp133 = -_tmp111 * _tmp127 - _tmp118 * _tmp95 - _tmp130 * _tmp132 -
                         _tmp35 * _tmp94 * (-_tmp45 * _tmp91 + _tmp76 * _tmp92);
  const Scalar _tmp134 = Scalar(1.0) / (_tmp133);
  const Scalar _tmp135 = _tmp50 + _tmp69;
  const Scalar _tmp136 = _tmp135 * _tmp53;
  const Scalar _tmp137 = Scalar(1.0) / (-_tmp136 + _tmp47 - _tmp72);
  const Scalar _tmp138 = Scalar(1.0) * _tmp137;
  const Scalar _tmp139 = _tmp105 * _tmp138;
  const Scalar _tmp140 = -_tmp131 * _tmp87 + _tmp135 * _tmp139;
  const Scalar _tmp141 = Scalar(1.0) * _tmp139 - Scalar(1.0) * _tmp140 * _tmp97;
  const Scalar _tmp142 = _tmp87 * _tmp89;
  const Scalar _tmp143 = _tmp135 * _tmp137;
  const Scalar _tmp144 = _tmp120 + _tmp122 * _tmp143 - _tmp123 * _tmp142;
  const Scalar _tmp145 = _tmp122 * _tmp138 - _tmp144 * _tmp97;
  const Scalar _tmp146 = Scalar(1.0) * _tmp116;
  const Scalar _tmp147 = _tmp145 * _tmp146;
  const Scalar _tmp148 = _tmp67 + _tmp77;
  const Scalar _tmp149 = _tmp116 * _tmp148;
  const Scalar _tmp150 = -_tmp149 * _tmp95 - Scalar(5.1796800000000003) * _tmp25 - _tmp30 * fv1;
  const Scalar _tmp151 = _tmp138 * _tmp53;
  const Scalar _tmp152 = _tmp136 * _tmp138 + Scalar(1.0);
  const Scalar _tmp153 = Scalar(1.0) * _tmp151 - Scalar(1.0) * _tmp152 * _tmp97;
  const Scalar _tmp154 = _tmp106 * _tmp143 - _tmp107 * _tmp142 - Scalar(1.0) * _tmp86;
  const Scalar _tmp155 = _tmp106 * _tmp138 - _tmp154 * _tmp97;
  const Scalar _tmp156 = _tmp146 * _tmp155;
  const Scalar _tmp157 = _tmp137 * _tmp54;
  const Scalar _tmp158 = -_tmp135 * _tmp157 - _tmp142 * _tmp55 + _tmp50;
  const Scalar _tmp159 = Scalar(5.1796800000000003) * _tmp11 + _tmp111 * _tmp149 + _tmp17 * fv1;
  const Scalar _tmp160 = _tmp135 * _tmp51;
  const Scalar _tmp161 = _tmp138 * _tmp160;
  const Scalar _tmp162 = -Scalar(1.0) * _tmp138 + Scalar(1.0) * _tmp161;
  const Scalar _tmp163 =
      _tmp111 * _tmp147 + _tmp130 * _tmp141 + _tmp150 * _tmp153 + _tmp156 * _tmp95 +
      _tmp159 * _tmp162 +
      Scalar(1.0) * _tmp35 * (-_tmp138 * _tmp54 - _tmp158 * _tmp97 + Scalar(1.0));
  const Scalar _tmp164 = std::asinh(_tmp134 * _tmp163);
  const Scalar _tmp165 = Scalar(1.0) * _tmp164;
  const Scalar _tmp166 = Scalar(1.0) * std::sinh(_tmp165);
  const Scalar _tmp167 = std::pow(_tmp133, Scalar(-2));
  const Scalar _tmp168 =
      std::pow(Scalar(std::pow(_tmp163, Scalar(2)) * _tmp167 + 1), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp169 = _tmp148 * fh1;
  const Scalar _tmp170 = std::pow(_tmp114, Scalar(Scalar(-3) / Scalar(2)));
  const Scalar _tmp171 = _tmp111 * _tmp95;
  const Scalar _tmp172 = _tmp170 * _tmp171;
  const Scalar _tmp173 = _tmp169 * _tmp172;
  const Scalar _tmp174 = _tmp112 * _tmp170;
  const Scalar _tmp175 = _tmp174 * fh1;
  const Scalar _tmp176 = Scalar(1.0) * _tmp145;
  const Scalar _tmp177 = _tmp170 * _tmp30;
  const Scalar _tmp178 = fh1 * (-_tmp112 * _tmp177 + _tmp129 + _tmp17 * _tmp172);
  const Scalar _tmp179 = -_tmp149 + _tmp169 * _tmp174;
  const Scalar _tmp180 = Scalar(1.0) * _tmp155;
  const Scalar _tmp181 = _tmp172 * fh1;
  const Scalar _tmp182 = _tmp110 * _tmp94;
  const Scalar _tmp183 = _tmp126 * _tmp94;
  const Scalar _tmp184 = _tmp127 - _tmp132 * _tmp178 - _tmp175 * _tmp183 - _tmp181 * _tmp182;
  const Scalar _tmp185 = _tmp163 * _tmp167;
  const Scalar _tmp186 =
      _tmp168 * (_tmp134 * (_tmp141 * _tmp178 - _tmp147 - _tmp153 * _tmp173 + _tmp162 * _tmp179 +
                            _tmp175 * _tmp176 + _tmp180 * _tmp181) -
                 _tmp184 * _tmp185);
  const Scalar _tmp187 = Scalar(9.6622558468725703) * _tmp164;
  const Scalar _tmp188 =
      -_tmp133 * _tmp187 - std::sqrt(Scalar(std::pow(Scalar(-_tmp56 + p_b(1, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp58 + p_b(0, 0)), Scalar(2))));
  const Scalar _tmp189 = Scalar(0.1034955) * _tmp134;
  const Scalar _tmp190 = _tmp188 * _tmp189;
  const Scalar _tmp191 = std::sinh(_tmp190);
  const Scalar _tmp192 = Scalar(0.1034955) * _tmp167;
  const Scalar _tmp193 = _tmp184 * _tmp192;
  const Scalar _tmp194 = Scalar(9.6622558468725703) * _tmp133;
  const Scalar _tmp195 = _tmp189 * p_b(2, 0) + std::cosh(_tmp165) - std::cosh(_tmp190);
  const Scalar _tmp196 = _tmp130 * _tmp131;
  const Scalar _tmp197 = _tmp116 * _tmp124;
  const Scalar _tmp198 = _tmp108 * _tmp116;
  const Scalar _tmp199 = _tmp111 * _tmp197 - _tmp196 * _tmp84 + _tmp198 * _tmp95 - _tmp35 * _tmp91;
  const Scalar _tmp200 = Scalar(1.0) / (_tmp199);
  const Scalar _tmp201 = _tmp152 * _tmp51;
  const Scalar _tmp202 = _tmp140 * _tmp51;
  const Scalar _tmp203 = _tmp138 * _tmp159;
  const Scalar _tmp204 = _tmp116 * _tmp51;
  const Scalar _tmp205 = _tmp144 * _tmp204;
  const Scalar _tmp206 = _tmp154 * _tmp204;
  const Scalar _tmp207 = _tmp111 * _tmp205 + _tmp130 * _tmp202 + _tmp150 * _tmp201 +
                         _tmp158 * _tmp35 * _tmp51 - _tmp160 * _tmp203 + _tmp206 * _tmp95;
  const Scalar _tmp208 = std::asinh(_tmp200 * _tmp207);
  const Scalar _tmp209 = Scalar(9.6622558468725703) * _tmp199;
  const Scalar _tmp210 =
      -_tmp208 * _tmp209 - std::sqrt(Scalar(std::pow(Scalar(-_tmp38 + p_c(0, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp42 + p_c(1, 0)), Scalar(2))));
  const Scalar _tmp211 = _tmp131 * _tmp178;
  const Scalar _tmp212 = _tmp108 * _tmp181 + _tmp124 * _tmp175 - _tmp197 - _tmp211 * _tmp84;
  const Scalar _tmp213 = std::pow(_tmp199, Scalar(-2));
  const Scalar _tmp214 = Scalar(0.1034955) * _tmp213;
  const Scalar _tmp215 = _tmp212 * _tmp214;
  const Scalar _tmp216 = _tmp207 * _tmp213;
  const Scalar _tmp217 = _tmp154 * _tmp51;
  const Scalar _tmp218 = _tmp144 * _tmp51;
  const Scalar _tmp219 = _tmp138 * _tmp179;
  const Scalar _tmp220 = _tmp200 * (-_tmp160 * _tmp219 - _tmp173 * _tmp201 + _tmp175 * _tmp218 +
                                    _tmp178 * _tmp202 + _tmp181 * _tmp217 - _tmp205) -
                         _tmp212 * _tmp216;
  const Scalar _tmp221 =
      std::pow(Scalar(std::pow(_tmp207, Scalar(2)) * _tmp213 + 1), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp222 = _tmp209 * _tmp221;
  const Scalar _tmp223 = Scalar(9.6622558468725703) * _tmp212;
  const Scalar _tmp224 = Scalar(0.1034955) * _tmp200;
  const Scalar _tmp225 = _tmp210 * _tmp224;
  const Scalar _tmp226 = std::sinh(_tmp225);
  const Scalar _tmp227 = Scalar(1.0) * _tmp208;
  const Scalar _tmp228 = Scalar(1.0) * _tmp221 * std::sinh(_tmp227);
  const Scalar _tmp229 = _tmp224 * p_c(2, 0) - std::cosh(_tmp225) + std::cosh(_tmp227);
  const Scalar _tmp230 = _tmp116 * _tmp125;
  const Scalar _tmp231 = _tmp109 * _tmp181 + _tmp125 * _tmp175 + _tmp211 - _tmp230;
  const Scalar _tmp232 = _tmp109 * _tmp116;
  const Scalar _tmp233 = _tmp111 * _tmp230 + _tmp196 + _tmp232 * _tmp95 + _tmp35 * _tmp92;
  const Scalar _tmp234 = std::pow(_tmp233, Scalar(-2));
  const Scalar _tmp235 = Scalar(0.1034955) * _tmp234;
  const Scalar _tmp236 = _tmp235 * p_d(2, 0);
  const Scalar _tmp237 = _tmp116 * _tmp137;
  const Scalar _tmp238 = _tmp106 * _tmp237;
  const Scalar _tmp239 = _tmp122 * _tmp237;
  const Scalar _tmp240 = -_tmp111 * _tmp239 - _tmp130 * _tmp139 - _tmp150 * _tmp151 +
                         _tmp157 * _tmp35 + _tmp203 - _tmp238 * _tmp95;
  const Scalar _tmp241 =
      std::pow(Scalar(_tmp234 * std::pow(_tmp240, Scalar(2)) + 1), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp242 = _tmp234 * _tmp240;
  const Scalar _tmp243 = Scalar(1.0) / (_tmp233);
  const Scalar _tmp244 = _tmp122 * _tmp137;
  const Scalar _tmp245 = _tmp106 * _tmp137;
  const Scalar _tmp246 =
      _tmp241 *
      (-_tmp231 * _tmp242 + _tmp243 * (-_tmp139 * _tmp178 + _tmp151 * _tmp173 - _tmp175 * _tmp244 -
                                       _tmp181 * _tmp245 + _tmp219 + _tmp239));
  const Scalar _tmp247 = std::asinh(_tmp240 * _tmp243);
  const Scalar _tmp248 = Scalar(1.0) * _tmp247;
  const Scalar _tmp249 = Scalar(1.0) * std::sinh(_tmp248);
  const Scalar _tmp250 = Scalar(9.6622558468725703) * _tmp247;
  const Scalar _tmp251 =
      -_tmp233 * _tmp250 - std::sqrt(Scalar(std::pow(Scalar(-_tmp70 + p_d(1, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp73 + p_d(0, 0)), Scalar(2))));
  const Scalar _tmp252 = Scalar(0.1034955) * _tmp243;
  const Scalar _tmp253 = _tmp251 * _tmp252;
  const Scalar _tmp254 = std::sinh(_tmp253);
  const Scalar _tmp255 = Scalar(9.6622558468725703) * _tmp233;
  const Scalar _tmp256 = _tmp235 * _tmp251;
  const Scalar _tmp257 = Scalar(9.6622558468725703) * _tmp252 * p_d(2, 0) +
                         Scalar(9.6622558468725703) * std::cosh(_tmp248) -
                         Scalar(9.6622558468725703) * std::cosh(_tmp253);
  const Scalar _tmp258 = _tmp113 * _tmp170;
  const Scalar _tmp259 = _tmp258 * fh1;
  const Scalar _tmp260 = fh1 * (-_tmp128 + _tmp17 * _tmp258 - _tmp171 * _tmp177);
  const Scalar _tmp261 = _tmp118 - _tmp132 * _tmp260 - _tmp181 * _tmp183 - _tmp182 * _tmp259;
  const Scalar _tmp262 = Scalar(9.6622558468725703) * _tmp261;
  const Scalar _tmp263 = _tmp192 * _tmp261;
  const Scalar _tmp264 = _tmp149 - _tmp169 * _tmp258;
  const Scalar _tmp265 =
      _tmp168 * (_tmp134 * (_tmp141 * _tmp260 + _tmp153 * _tmp264 - _tmp156 + _tmp162 * _tmp173 +
                            _tmp176 * _tmp181 + _tmp180 * _tmp259) -
                 _tmp185 * _tmp261);
  const Scalar _tmp266 = _tmp131 * _tmp260;
  const Scalar _tmp267 = _tmp108 * _tmp259 + _tmp124 * _tmp181 - _tmp198 - _tmp266 * _tmp84;
  const Scalar _tmp268 = Scalar(9.6622558468725703) * _tmp267;
  const Scalar _tmp269 = _tmp214 * _tmp267;
  const Scalar _tmp270 = _tmp200 * (-_tmp161 * _tmp173 + _tmp181 * _tmp218 + _tmp201 * _tmp264 +
                                    _tmp202 * _tmp260 - _tmp206 + _tmp217 * _tmp259) -
                         _tmp216 * _tmp267;
  const Scalar _tmp271 = _tmp109 * _tmp259 + _tmp125 * _tmp181 - _tmp232 + _tmp266;
  const Scalar _tmp272 =
      _tmp241 *
      (-_tmp242 * _tmp271 + _tmp243 * (_tmp138 * _tmp173 - _tmp139 * _tmp260 - _tmp151 * _tmp264 -
                                       _tmp181 * _tmp244 + _tmp238 - _tmp245 * _tmp259));

  // Output terms (1)
  Eigen::Matrix<Scalar, 4, 3> _res;

  _res(0, 0) = -_tmp19 * _tmp34;
  _res(1, 0) =
      -Scalar(9.6622558468725703) * _tmp184 * _tmp195 -
      _tmp194 *
          (_tmp166 * _tmp186 -
           _tmp191 * (-_tmp188 * _tmp193 + _tmp189 * (-_tmp184 * _tmp187 - _tmp186 * _tmp194)) -
           _tmp193 * p_b(2, 0));
  _res(2, 0) =
      -_tmp209 *
          (-_tmp215 * p_c(2, 0) + _tmp220 * _tmp228 -
           _tmp226 * (-_tmp210 * _tmp215 + _tmp224 * (-_tmp208 * _tmp223 - _tmp220 * _tmp222))) -
      _tmp223 * _tmp229;
  _res(3, 0) =
      -_tmp231 * _tmp257 -
      _tmp255 *
          (-_tmp231 * _tmp236 + _tmp246 * _tmp249 -
           _tmp254 * (-_tmp231 * _tmp256 + _tmp252 * (-_tmp231 * _tmp250 - _tmp246 * _tmp255)));
  _res(0, 1) = -_tmp32 * _tmp34;
  _res(1, 1) =
      -_tmp194 *
          (_tmp166 * _tmp265 -
           _tmp191 * (-_tmp188 * _tmp263 + _tmp189 * (-_tmp164 * _tmp262 - _tmp194 * _tmp265)) -
           _tmp263 * p_b(2, 0)) -
      _tmp195 * _tmp262;
  _res(2, 1) =
      -_tmp209 *
          (-_tmp226 * (-_tmp210 * _tmp269 + _tmp224 * (-_tmp208 * _tmp268 - _tmp222 * _tmp270)) +
           _tmp228 * _tmp270 - _tmp269 * p_c(2, 0)) -
      _tmp229 * _tmp268;
  _res(3, 1) =
      -_tmp255 *
          (-_tmp236 * _tmp271 + _tmp249 * _tmp272 -
           _tmp254 * (_tmp252 * (-_tmp250 * _tmp271 - _tmp255 * _tmp272) - _tmp256 * _tmp271)) -
      _tmp257 * _tmp271;
  _res(0, 2) = Scalar(-1.0);
  _res(1, 2) = 0;
  _res(2, 2) = 0;
  _res(3, 2) = 0;

  return _res;
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym

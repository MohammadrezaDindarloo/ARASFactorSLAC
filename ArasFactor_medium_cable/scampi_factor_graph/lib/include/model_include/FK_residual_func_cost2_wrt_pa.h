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
 * Symbolic function: FK_residual_func_cost2_wrt_pa
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
Eigen::Matrix<Scalar, 4, 3> FkResidualFuncCost2WrtPa(
    const Scalar fh1, const Scalar fv1, const sym::Rot3<Scalar>& DeltaRot,
    const sym::Pose3<Scalar>& TransformationMatrix, const Eigen::Matrix<Scalar, 4, 1>& encoder,
    const Eigen::Matrix<Scalar, 3, 1>& p_a, const Eigen::Matrix<Scalar, 3, 1>& p_b,
    const Eigen::Matrix<Scalar, 3, 1>& p_c, const Eigen::Matrix<Scalar, 3, 1>& p_d,
    const Scalar epsilon) {
  // Total ops: 767

  // Unused inputs
  (void)encoder;
  (void)epsilon;

  // Input arrays
  const Eigen::Matrix<Scalar, 4, 1>& _DeltaRot = DeltaRot.Data();
  const Eigen::Matrix<Scalar, 7, 1>& _TransformationMatrix = TransformationMatrix.Data();

  // Intermediate terms (264)
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
  const Scalar _tmp13 = -2 * std::pow(_tmp1, Scalar(2));
  const Scalar _tmp14 = 1 - 2 * std::pow(_tmp3, Scalar(2));
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
  const Scalar _tmp28 = Scalar(0.20999999999999999) * _tmp14 + Scalar(0.20999999999999999) * _tmp27;
  const Scalar _tmp29 = _tmp26 - _tmp28;
  const Scalar _tmp30 = _tmp22 + _tmp29;
  const Scalar _tmp31 = _TransformationMatrix[5] + _tmp30;
  const Scalar _tmp32 = -_tmp31 + p_a(1, 0);
  const Scalar _tmp33 = std::pow(_tmp19, Scalar(2)) + std::pow(_tmp32, Scalar(2));
  const Scalar _tmp34 = std::sqrt(_tmp33);
  const Scalar _tmp35 =
      Scalar(1.0) *
      std::cosh(Scalar(0.1034955) * _tmp20 *
                (-_tmp34 - Scalar(9.6622558468725703) * fh1 * std::asinh(_tmp20 * fv1))) /
      _tmp34;
  const Scalar _tmp36 = Scalar(0.20999999999999999) * _tmp23 + Scalar(0.20999999999999999) * _tmp24;
  const Scalar _tmp37 = -_tmp36;
  const Scalar _tmp38 = -Scalar(0.010999999999999999) * _tmp13 -
                        Scalar(0.010999999999999999) * _tmp27 + Scalar(-0.010999999999999999);
  const Scalar _tmp39 = -Scalar(0.20999999999999999) * _tmp10 + Scalar(0.20999999999999999) * _tmp9;
  const Scalar _tmp40 = _tmp38 - _tmp39;
  const Scalar _tmp41 = _tmp37 + _tmp40;
  const Scalar _tmp42 = -_TransformationMatrix[6] - _tmp41 + p_a(2, 0);
  const Scalar _tmp43 =
      std::pow(Scalar(_tmp33 + std::pow(_tmp42, Scalar(2))), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp44 = _tmp18 - p_a(0, 0);
  const Scalar _tmp45 = _tmp31 - p_a(1, 0);
  const Scalar _tmp46 = _tmp44 * _tmp45;
  const Scalar _tmp47 = std::pow(_tmp44, Scalar(2));
  const Scalar _tmp48 = std::pow(_tmp45, Scalar(2));
  const Scalar _tmp49 = _tmp47 + _tmp48;
  const Scalar _tmp50 = std::pow(_tmp49, Scalar(Scalar(-3) / Scalar(2)));
  const Scalar _tmp51 = _tmp50 * fh1;
  const Scalar _tmp52 = _tmp46 * _tmp51;
  const Scalar _tmp53 = _tmp12 + _tmp15;
  const Scalar _tmp54 = _tmp53 + _tmp6;
  const Scalar _tmp55 = _TransformationMatrix[4] + _tmp54;
  const Scalar _tmp56 = _tmp55 - p_c(0, 0);
  const Scalar _tmp57 = _tmp26 + _tmp28;
  const Scalar _tmp58 = _tmp21 + _tmp57;
  const Scalar _tmp59 = _TransformationMatrix[5] + _tmp58;
  const Scalar _tmp60 = _tmp59 - p_c(1, 0);
  const Scalar _tmp61 = std::pow(Scalar(std::pow(_tmp56, Scalar(2)) + std::pow(_tmp60, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp62 = _tmp56 * _tmp61;
  const Scalar _tmp63 = _tmp21 + _tmp29;
  const Scalar _tmp64 = _TransformationMatrix[5] + _tmp63;
  const Scalar _tmp65 = _tmp64 - p_b(1, 0);
  const Scalar _tmp66 = _tmp53 + _tmp7;
  const Scalar _tmp67 = _TransformationMatrix[4] + _tmp66;
  const Scalar _tmp68 = _tmp67 - p_b(0, 0);
  const Scalar _tmp69 = Scalar(1.0) / (_tmp68);
  const Scalar _tmp70 = _tmp65 * _tmp69;
  const Scalar _tmp71 = _tmp60 * _tmp61;
  const Scalar _tmp72 = Scalar(1.0) / (_tmp62 * _tmp70 - _tmp71);
  const Scalar _tmp73 = _tmp38 + _tmp39;
  const Scalar _tmp74 = _tmp36 + _tmp73;
  const Scalar _tmp75 = _tmp37 + _tmp73;
  const Scalar _tmp76 = _tmp72 * (-_tmp62 * _tmp74 + _tmp62 * _tmp75);
  const Scalar _tmp77 = Scalar(1.0) * _tmp66;
  const Scalar _tmp78 = -_tmp54 + _tmp77;
  const Scalar _tmp79 = Scalar(1.0) * _tmp63;
  const Scalar _tmp80 = -_tmp79;
  const Scalar _tmp81 = Scalar(1.0) / (_tmp58 + _tmp80);
  const Scalar _tmp82 = Scalar(1.0) * _tmp81;
  const Scalar _tmp83 = _tmp70 * _tmp75;
  const Scalar _tmp84 = _tmp72 * (-_tmp62 * _tmp83 + _tmp71 * _tmp74);
  const Scalar _tmp85 = -Scalar(1.0) * _tmp76 + _tmp78 * _tmp82 * _tmp84;
  const Scalar _tmp86 = _tmp36 + _tmp40;
  const Scalar _tmp87 = _tmp22 + _tmp57;
  const Scalar _tmp88 = _TransformationMatrix[5] + _tmp87;
  const Scalar _tmp89 = _tmp88 - p_d(1, 0);
  const Scalar _tmp90 = _tmp16 + _tmp6;
  const Scalar _tmp91 = _TransformationMatrix[4] + _tmp90;
  const Scalar _tmp92 = _tmp91 - p_d(0, 0);
  const Scalar _tmp93 = std::pow(Scalar(std::pow(_tmp89, Scalar(2)) + std::pow(_tmp92, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp94 = _tmp92 * _tmp93;
  const Scalar _tmp95 = _tmp89 * _tmp93;
  const Scalar _tmp96 = _tmp70 * _tmp94 - _tmp95;
  const Scalar _tmp97 = -_tmp83 * _tmp94 - _tmp84 * _tmp96 + _tmp86 * _tmp95;
  const Scalar _tmp98 = _tmp78 * _tmp81;
  const Scalar _tmp99 = _tmp75 * _tmp94 - _tmp76 * _tmp96 - _tmp86 * _tmp94 - _tmp97 * _tmp98;
  const Scalar _tmp100 = Scalar(1.0) / (_tmp99);
  const Scalar _tmp101 =
      std::sqrt(Scalar(std::pow(_tmp65, Scalar(2)) + std::pow(_tmp68, Scalar(2))));
  const Scalar _tmp102 = Scalar(1.0) / (_tmp101);
  const Scalar _tmp103 = _tmp101 * _tmp69;
  const Scalar _tmp104 = _tmp103 * (-_tmp102 * _tmp63 * _tmp68 + _tmp102 * _tmp65 * _tmp66);
  const Scalar _tmp105 = _tmp72 * (_tmp104 * _tmp62 - _tmp54 * _tmp71 + _tmp58 * _tmp62);
  const Scalar _tmp106 = _tmp104 * _tmp94 - _tmp105 * _tmp96 + _tmp87 * _tmp94 - _tmp90 * _tmp95;
  const Scalar _tmp107 = _tmp100 * _tmp106;
  const Scalar _tmp108 = Scalar(1.0) / (_tmp106);
  const Scalar _tmp109 = _tmp108 * _tmp99;
  const Scalar _tmp110 = _tmp109 * (-Scalar(1.0) * _tmp105 - _tmp107 * _tmp85);
  const Scalar _tmp111 = _tmp110 + _tmp85;
  const Scalar _tmp112 = _tmp100 * _tmp96;
  const Scalar _tmp113 = _tmp72 * (-_tmp111 * _tmp112 + Scalar(1.0));
  const Scalar _tmp114 = _tmp100 * _tmp111;
  const Scalar _tmp115 = _tmp113 * _tmp62 + _tmp114 * _tmp94;
  const Scalar _tmp116 = _tmp103 * _tmp115;
  const Scalar _tmp117 = _tmp47 * _tmp51;
  const Scalar _tmp118 = _tmp70 * _tmp84 + _tmp83;
  const Scalar _tmp119 = -_tmp118 * _tmp98 + _tmp70 * _tmp76 - _tmp75;
  const Scalar _tmp120 = _tmp109 * (-_tmp104 + _tmp105 * _tmp70 - _tmp107 * _tmp119);
  const Scalar _tmp121 = _tmp119 + _tmp120;
  const Scalar _tmp122 = _tmp72 * (-_tmp112 * _tmp121 - _tmp70);
  const Scalar _tmp123 = _tmp100 * _tmp121;
  const Scalar _tmp124 = _tmp122 * _tmp62 + _tmp123 * _tmp94 + Scalar(1.0);
  const Scalar _tmp125 = _tmp103 * _tmp124;
  const Scalar _tmp126 = std::pow(_tmp49, Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp127 = _tmp126 * fh1;
  const Scalar _tmp128 = _tmp103 * _tmp127;
  const Scalar _tmp129 = _tmp124 * _tmp128;
  const Scalar _tmp130 = _tmp126 * _tmp30;
  const Scalar _tmp131 = _tmp30 * _tmp50;
  const Scalar _tmp132 = _tmp17 * _tmp50;
  const Scalar _tmp133 = fh1 * (_tmp130 - _tmp131 * _tmp47 + _tmp132 * _tmp46);
  const Scalar _tmp134 = Scalar(1.0) * _tmp108;
  const Scalar _tmp135 = _tmp72 * _tmp96;
  const Scalar _tmp136 = _tmp103 * (-_tmp134 * _tmp135 * _tmp62 + _tmp134 * _tmp94);
  const Scalar _tmp137 = -_tmp116 * _tmp52 - _tmp117 * _tmp125 + _tmp129 - _tmp133 * _tmp136;
  const Scalar _tmp138 = Scalar(43.164000000000001) - fv1;
  const Scalar _tmp139 = _tmp77 + _tmp79 * _tmp98;
  const Scalar _tmp140 = 0;
  const Scalar _tmp141 = _tmp112 * _tmp140 * _tmp72;
  const Scalar _tmp142 = _tmp100 * _tmp140;
  const Scalar _tmp143 = _tmp115 * _tmp128;
  const Scalar _tmp144 = _tmp126 * _tmp17;
  const Scalar _tmp145 = fh1 * (-_tmp130 * _tmp44 + _tmp144 * _tmp45);
  const Scalar _tmp146 = -_tmp103 * _tmp138 * (-_tmp141 * _tmp62 + _tmp142 * _tmp94) -
                         _tmp129 * _tmp44 - _tmp136 * _tmp145 - _tmp143 * _tmp45;
  const Scalar _tmp147 = Scalar(1.0) / (_tmp146);
  const Scalar _tmp148 = _tmp80 + _tmp87;
  const Scalar _tmp149 = _tmp148 * _tmp98;
  const Scalar _tmp150 = Scalar(1.0) / (-_tmp149 + _tmp77 - _tmp90);
  const Scalar _tmp151 = Scalar(1.0) * _tmp150;
  const Scalar _tmp152 = _tmp109 * _tmp151;
  const Scalar _tmp153 = -_tmp134 * _tmp97 + _tmp148 * _tmp152;
  const Scalar _tmp154 = Scalar(1.0) * _tmp152 - Scalar(1.0) * _tmp153 * _tmp82;
  const Scalar _tmp155 = _tmp100 * _tmp97;
  const Scalar _tmp156 = _tmp148 * _tmp150;
  const Scalar _tmp157 = _tmp118 + _tmp120 * _tmp156 - _tmp121 * _tmp155;
  const Scalar _tmp158 = _tmp120 * _tmp151 - _tmp157 * _tmp82;
  const Scalar _tmp159 = Scalar(1.0) * _tmp127;
  const Scalar _tmp160 = _tmp158 * _tmp159;
  const Scalar _tmp161 = _tmp127 * _tmp41;
  const Scalar _tmp162 = -_tmp161 * _tmp45 - Scalar(5.1796800000000003) * _tmp25 - _tmp30 * fv1;
  const Scalar _tmp163 = _tmp151 * _tmp98;
  const Scalar _tmp164 = _tmp149 * _tmp151 + Scalar(1.0);
  const Scalar _tmp165 = Scalar(1.0) * _tmp163 - Scalar(1.0) * _tmp164 * _tmp82;
  const Scalar _tmp166 = _tmp110 * _tmp156 - _tmp111 * _tmp155 - Scalar(1.0) * _tmp84;
  const Scalar _tmp167 = _tmp110 * _tmp151 - _tmp166 * _tmp82;
  const Scalar _tmp168 = _tmp159 * _tmp167;
  const Scalar _tmp169 = _tmp139 * _tmp150;
  const Scalar _tmp170 = -_tmp140 * _tmp155 - _tmp148 * _tmp169 + _tmp80;
  const Scalar _tmp171 = Scalar(5.1796800000000003) * _tmp11 + _tmp161 * _tmp44 + _tmp17 * fv1;
  const Scalar _tmp172 = _tmp148 * _tmp81;
  const Scalar _tmp173 = _tmp151 * _tmp172;
  const Scalar _tmp174 = -Scalar(1.0) * _tmp151 + Scalar(1.0) * _tmp173;
  const Scalar _tmp175 =
      Scalar(1.0) * _tmp138 * (-_tmp139 * _tmp151 - _tmp170 * _tmp82 + Scalar(1.0)) +
      _tmp145 * _tmp154 + _tmp160 * _tmp44 + _tmp162 * _tmp165 + _tmp168 * _tmp45 +
      _tmp171 * _tmp174;
  const Scalar _tmp176 = std::asinh(_tmp147 * _tmp175);
  const Scalar _tmp177 = Scalar(9.6622558468725703) * _tmp176;
  const Scalar _tmp178 =
      -_tmp146 * _tmp177 - std::sqrt(Scalar(std::pow(Scalar(-_tmp64 + p_b(1, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp67 + p_b(0, 0)), Scalar(2))));
  const Scalar _tmp179 = Scalar(0.1034955) * _tmp147;
  const Scalar _tmp180 = _tmp178 * _tmp179;
  const Scalar _tmp181 = Scalar(1.0) * _tmp176;
  const Scalar _tmp182 = -Scalar(9.6622558468725703) * std::sinh(_tmp180) -
                         Scalar(9.6622558468725703) * std::sinh(_tmp181);
  const Scalar _tmp183 = Scalar(1.0) * std::cosh(_tmp181);
  const Scalar _tmp184 = std::pow(_tmp146, Scalar(-2));
  const Scalar _tmp185 =
      std::pow(Scalar(std::pow(_tmp175, Scalar(2)) * _tmp184 + 1), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp186 = _tmp41 * _tmp51;
  const Scalar _tmp187 = _tmp186 * _tmp46;
  const Scalar _tmp188 = Scalar(1.0) * _tmp158;
  const Scalar _tmp189 = -_tmp161 + _tmp186 * _tmp47;
  const Scalar _tmp190 = Scalar(1.0) * _tmp167;
  const Scalar _tmp191 = _tmp175 * _tmp184;
  const Scalar _tmp192 =
      _tmp185 *
      (-_tmp137 * _tmp191 + _tmp147 * (_tmp117 * _tmp188 + _tmp133 * _tmp154 - _tmp160 -
                                       _tmp165 * _tmp187 + _tmp174 * _tmp189 + _tmp190 * _tmp52));
  const Scalar _tmp193 = std::cosh(_tmp180);
  const Scalar _tmp194 = Scalar(0.1034955) * _tmp178 * _tmp184;
  const Scalar _tmp195 = Scalar(9.6622558468725703) * _tmp146;
  const Scalar _tmp196 = _tmp122 * _tmp127;
  const Scalar _tmp197 = _tmp133 * _tmp134;
  const Scalar _tmp198 = _tmp113 * _tmp52 + _tmp117 * _tmp122 - _tmp135 * _tmp197 - _tmp196;
  const Scalar _tmp199 = _tmp134 * _tmp145;
  const Scalar _tmp200 = _tmp113 * _tmp127;
  const Scalar _tmp201 =
      -_tmp135 * _tmp199 - _tmp138 * _tmp141 + _tmp196 * _tmp44 + _tmp200 * _tmp45;
  const Scalar _tmp202 = Scalar(1.0) / (_tmp201);
  const Scalar _tmp203 = _tmp164 * _tmp81;
  const Scalar _tmp204 = _tmp153 * _tmp81;
  const Scalar _tmp205 = _tmp151 * _tmp171;
  const Scalar _tmp206 = _tmp127 * _tmp81;
  const Scalar _tmp207 = _tmp157 * _tmp206;
  const Scalar _tmp208 = _tmp166 * _tmp206;
  const Scalar _tmp209 = _tmp138 * _tmp170 * _tmp81 + _tmp145 * _tmp204 + _tmp162 * _tmp203 -
                         _tmp172 * _tmp205 + _tmp207 * _tmp44 + _tmp208 * _tmp45;
  const Scalar _tmp210 = std::asinh(_tmp202 * _tmp209);
  const Scalar _tmp211 = Scalar(1.0) * _tmp210;
  const Scalar _tmp212 = Scalar(9.6622558468725703) * _tmp210;
  const Scalar _tmp213 =
      -_tmp201 * _tmp212 - std::sqrt(Scalar(std::pow(Scalar(-_tmp55 + p_c(0, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp59 + p_c(1, 0)), Scalar(2))));
  const Scalar _tmp214 = Scalar(0.1034955) * _tmp202;
  const Scalar _tmp215 = _tmp213 * _tmp214;
  const Scalar _tmp216 = -Scalar(9.6622558468725703) * std::sinh(_tmp211) -
                         Scalar(9.6622558468725703) * std::sinh(_tmp215);
  const Scalar _tmp217 = std::pow(_tmp201, Scalar(-2));
  const Scalar _tmp218 = Scalar(0.1034955) * _tmp213 * _tmp217;
  const Scalar _tmp219 = Scalar(9.6622558468725703) * _tmp201;
  const Scalar _tmp220 = _tmp209 * _tmp217;
  const Scalar _tmp221 = _tmp166 * _tmp81;
  const Scalar _tmp222 = _tmp157 * _tmp81;
  const Scalar _tmp223 = _tmp151 * _tmp189;
  const Scalar _tmp224 =
      std::pow(Scalar(std::pow(_tmp209, Scalar(2)) * _tmp217 + 1), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp225 =
      _tmp224 *
      (-_tmp198 * _tmp220 + _tmp202 * (_tmp117 * _tmp222 + _tmp133 * _tmp204 - _tmp172 * _tmp223 -
                                       _tmp187 * _tmp203 - _tmp207 + _tmp221 * _tmp52));
  const Scalar _tmp226 = std::cosh(_tmp215);
  const Scalar _tmp227 = Scalar(1.0) * std::cosh(_tmp211);
  const Scalar _tmp228 = _tmp123 * _tmp127;
  const Scalar _tmp229 = _tmp114 * _tmp52 + _tmp117 * _tmp123 + _tmp197 - _tmp228;
  const Scalar _tmp230 = _tmp114 * _tmp127;
  const Scalar _tmp231 = _tmp138 * _tmp142 + _tmp199 + _tmp228 * _tmp44 + _tmp230 * _tmp45;
  const Scalar _tmp232 = Scalar(1.0) / (_tmp231);
  const Scalar _tmp233 = _tmp127 * _tmp150;
  const Scalar _tmp234 = _tmp110 * _tmp233;
  const Scalar _tmp235 = _tmp120 * _tmp233;
  const Scalar _tmp236 = _tmp138 * _tmp169 - _tmp145 * _tmp152 - _tmp162 * _tmp163 + _tmp205 -
                         _tmp234 * _tmp45 - _tmp235 * _tmp44;
  const Scalar _tmp237 = std::asinh(_tmp232 * _tmp236);
  const Scalar _tmp238 = Scalar(9.6622558468725703) * _tmp231;
  const Scalar _tmp239 =
      -_tmp237 * _tmp238 - std::sqrt(Scalar(std::pow(Scalar(-_tmp88 + p_d(1, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp91 + p_d(0, 0)), Scalar(2))));
  const Scalar _tmp240 = Scalar(0.1034955) * _tmp232;
  const Scalar _tmp241 = _tmp239 * _tmp240;
  const Scalar _tmp242 = Scalar(1.0) * _tmp237;
  const Scalar _tmp243 = -Scalar(9.6622558468725703) * std::sinh(_tmp241) -
                         Scalar(9.6622558468725703) * std::sinh(_tmp242);
  const Scalar _tmp244 = std::pow(_tmp231, Scalar(-2));
  const Scalar _tmp245 =
      std::pow(Scalar(std::pow(_tmp236, Scalar(2)) * _tmp244 + 1), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp246 = _tmp236 * _tmp244;
  const Scalar _tmp247 = _tmp120 * _tmp150;
  const Scalar _tmp248 = _tmp110 * _tmp150;
  const Scalar _tmp249 =
      _tmp245 *
      (-_tmp229 * _tmp246 + _tmp232 * (-_tmp117 * _tmp247 - _tmp133 * _tmp152 + _tmp163 * _tmp187 +
                                       _tmp223 + _tmp235 - _tmp248 * _tmp52));
  const Scalar _tmp250 = Scalar(1.0) * std::cosh(_tmp242);
  const Scalar _tmp251 = std::cosh(_tmp241);
  const Scalar _tmp252 = Scalar(9.6622558468725703) * _tmp237;
  const Scalar _tmp253 = Scalar(0.1034955) * _tmp239 * _tmp244;
  const Scalar _tmp254 = _tmp161 - _tmp186 * _tmp48;
  const Scalar _tmp255 = _tmp48 * _tmp51;
  const Scalar _tmp256 = fh1 * (-_tmp131 * _tmp46 + _tmp132 * _tmp48 - _tmp144);
  const Scalar _tmp257 = -_tmp116 * _tmp255 - _tmp125 * _tmp52 - _tmp136 * _tmp256 + _tmp143;
  const Scalar _tmp258 =
      _tmp185 * (_tmp147 * (_tmp154 * _tmp256 + _tmp165 * _tmp254 - _tmp168 + _tmp174 * _tmp187 +
                            _tmp188 * _tmp52 + _tmp190 * _tmp255) -
                 _tmp191 * _tmp257);
  const Scalar _tmp259 = _tmp134 * _tmp256;
  const Scalar _tmp260 = _tmp113 * _tmp255 + _tmp122 * _tmp52 - _tmp135 * _tmp259 - _tmp200;
  const Scalar _tmp261 =
      _tmp224 * (_tmp202 * (-_tmp173 * _tmp187 + _tmp203 * _tmp254 + _tmp204 * _tmp256 - _tmp208 +
                            _tmp221 * _tmp255 + _tmp222 * _tmp52) -
                 _tmp220 * _tmp260);
  const Scalar _tmp262 = _tmp114 * _tmp255 + _tmp123 * _tmp52 - _tmp230 + _tmp259;
  const Scalar _tmp263 =
      _tmp245 * (_tmp232 * (_tmp151 * _tmp187 - _tmp152 * _tmp256 - _tmp163 * _tmp254 + _tmp234 -
                            _tmp247 * _tmp52 - _tmp248 * _tmp255) -
                 _tmp246 * _tmp262);

  // Output terms (1)
  Eigen::Matrix<Scalar, 4, 3> _res;

  _res(0, 0) = _tmp19 * _tmp35 - _tmp19 * _tmp43;
  _res(1, 0) =
      _tmp137 * _tmp182 +
      _tmp195 *
          (-_tmp183 * _tmp192 -
           _tmp193 * (-_tmp137 * _tmp194 + _tmp179 * (-_tmp137 * _tmp177 - _tmp192 * _tmp195)));
  _res(2, 0) =
      _tmp198 * _tmp216 +
      _tmp219 *
          (-_tmp225 * _tmp227 -
           _tmp226 * (-_tmp198 * _tmp218 + _tmp214 * (-_tmp198 * _tmp212 - _tmp219 * _tmp225)));
  _res(3, 0) =
      _tmp229 * _tmp243 +
      _tmp238 *
          (-_tmp249 * _tmp250 -
           _tmp251 * (-_tmp229 * _tmp253 + _tmp240 * (-_tmp229 * _tmp252 - _tmp238 * _tmp249)));
  _res(0, 1) = _tmp32 * _tmp35 - _tmp32 * _tmp43;
  _res(1, 1) =
      _tmp182 * _tmp257 + _tmp195 * (-_tmp183 * _tmp258 -
                                     _tmp193 * (_tmp179 * (-_tmp177 * _tmp257 - _tmp195 * _tmp258) -
                                                _tmp194 * _tmp257));
  _res(2, 1) =
      _tmp216 * _tmp260 +
      _tmp219 *
          (-_tmp226 * (_tmp214 * (-_tmp212 * _tmp260 - _tmp219 * _tmp261) - _tmp218 * _tmp260) -
           _tmp227 * _tmp261);
  _res(3, 1) =
      _tmp238 *
          (-_tmp250 * _tmp263 -
           _tmp251 * (_tmp240 * (-_tmp238 * _tmp263 - _tmp252 * _tmp262) - _tmp253 * _tmp262)) +
      _tmp243 * _tmp262;
  _res(0, 2) = -_tmp42 * _tmp43;
  _res(1, 2) = 0;
  _res(2, 2) = 0;
  _res(3, 2) = 0;

  return _res;
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym

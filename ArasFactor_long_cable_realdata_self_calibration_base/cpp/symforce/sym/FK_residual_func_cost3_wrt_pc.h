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
 * Symbolic function: FK_residual_func_cost3_wrt_pc
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
Eigen::Matrix<Scalar, 4, 3> FkResidualFuncCost3WrtPc(
    const Scalar fh1, const Scalar fv1, const sym::Rot3<Scalar>& DeltaRot,
    const sym::Pose3<Scalar>& TransformationMatrix, const Eigen::Matrix<Scalar, 4, 1>& encoder,
    const Eigen::Matrix<Scalar, 4, 1>& offset, const Eigen::Matrix<Scalar, 3, 1>& p_a,
    const Eigen::Matrix<Scalar, 3, 1>& p_b, const Eigen::Matrix<Scalar, 3, 1>& p_c,
    const Eigen::Matrix<Scalar, 3, 1>& p_d, const Scalar epsilon) {
  // Total ops: 1136

  // Unused inputs
  (void)encoder;
  (void)offset;
  (void)epsilon;

  // Input arrays
  const Eigen::Matrix<Scalar, 4, 1>& _DeltaRot = DeltaRot.Data();
  const Eigen::Matrix<Scalar, 7, 1>& _TransformationMatrix = TransformationMatrix.Data();

  // Intermediate terms (353)
  const Scalar _tmp0 =
      _DeltaRot[0] * _TransformationMatrix[3] - _DeltaRot[1] * _TransformationMatrix[2] +
      _DeltaRot[2] * _TransformationMatrix[1] + _DeltaRot[3] * _TransformationMatrix[0];
  const Scalar _tmp1 =
      _DeltaRot[0] * _TransformationMatrix[2] + _DeltaRot[1] * _TransformationMatrix[3] -
      _DeltaRot[2] * _TransformationMatrix[0] + _DeltaRot[3] * _TransformationMatrix[1];
  const Scalar _tmp2 = 2 * _tmp1;
  const Scalar _tmp3 = _tmp0 * _tmp2;
  const Scalar _tmp4 =
      -_DeltaRot[0] * _TransformationMatrix[1] + _DeltaRot[1] * _TransformationMatrix[0] +
      _DeltaRot[2] * _TransformationMatrix[3] + _DeltaRot[3] * _TransformationMatrix[2];
  const Scalar _tmp5 =
      -_DeltaRot[0] * _TransformationMatrix[0] - _DeltaRot[1] * _TransformationMatrix[1] -
      _DeltaRot[2] * _TransformationMatrix[2] + _DeltaRot[3] * _TransformationMatrix[3];
  const Scalar _tmp6 = 2 * _tmp5;
  const Scalar _tmp7 = _tmp4 * _tmp6;
  const Scalar _tmp8 = Scalar(0.20999999999999999) * _tmp3 - Scalar(0.20999999999999999) * _tmp7;
  const Scalar _tmp9 = -2 * std::pow(_tmp1, Scalar(2));
  const Scalar _tmp10 = 1 - 2 * std::pow(_tmp4, Scalar(2));
  const Scalar _tmp11 = Scalar(0.20999999999999999) * _tmp10 + Scalar(0.20999999999999999) * _tmp9;
  const Scalar _tmp12 = 2 * _tmp0 * _tmp4;
  const Scalar _tmp13 = _tmp2 * _tmp5;
  const Scalar _tmp14 = _tmp12 + _tmp13;
  const Scalar _tmp15 = -Scalar(0.010999999999999999) * _tmp14;
  const Scalar _tmp16 = _tmp11 + _tmp15;
  const Scalar _tmp17 = _tmp16 + _tmp8;
  const Scalar _tmp18 = _TransformationMatrix[4] + _tmp17;
  const Scalar _tmp19 = _tmp18 - p_c(0, 0);
  const Scalar _tmp20 = std::pow(_tmp19, Scalar(2));
  const Scalar _tmp21 = -2 * std::pow(_tmp0, Scalar(2));
  const Scalar _tmp22 = Scalar(0.20999999999999999) * _tmp10 + Scalar(0.20999999999999999) * _tmp21;
  const Scalar _tmp23 = Scalar(0.20999999999999999) * _tmp3 + Scalar(0.20999999999999999) * _tmp7;
  const Scalar _tmp24 = _tmp2 * _tmp4;
  const Scalar _tmp25 = _tmp0 * _tmp6;
  const Scalar _tmp26 = _tmp24 - _tmp25;
  const Scalar _tmp27 = -Scalar(0.010999999999999999) * _tmp26;
  const Scalar _tmp28 = _tmp23 + _tmp27;
  const Scalar _tmp29 = _tmp22 + _tmp28;
  const Scalar _tmp30 = _TransformationMatrix[5] + _tmp29;
  const Scalar _tmp31 = _tmp30 - p_c(1, 0);
  const Scalar _tmp32 = std::pow(_tmp31, Scalar(2));
  const Scalar _tmp33 = _tmp20 + _tmp32;
  const Scalar _tmp34 = std::pow(_tmp33, Scalar(Scalar(-3) / Scalar(2)));
  const Scalar _tmp35 = _tmp20 * _tmp34;
  const Scalar _tmp36 = -_tmp11 + _tmp15;
  const Scalar _tmp37 = _tmp36 + _tmp8;
  const Scalar _tmp38 = -_tmp23 + _tmp27;
  const Scalar _tmp39 = _tmp22 + _tmp38;
  const Scalar _tmp40 = _TransformationMatrix[5] + _tmp39;
  const Scalar _tmp41 = _tmp40 - p_d(1, 0);
  const Scalar _tmp42 = _TransformationMatrix[4] + _tmp37;
  const Scalar _tmp43 = _tmp42 - p_d(0, 0);
  const Scalar _tmp44 = std::pow(Scalar(std::pow(_tmp41, Scalar(2)) + std::pow(_tmp43, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp45 = _tmp41 * _tmp44;
  const Scalar _tmp46 = _tmp43 * _tmp44;
  const Scalar _tmp47 = std::pow(_tmp33, Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp48 = _tmp31 * _tmp47;
  const Scalar _tmp49 = _tmp29 * _tmp47;
  const Scalar _tmp50 = -_tmp8;
  const Scalar _tmp51 = _tmp16 + _tmp50;
  const Scalar _tmp52 = -_tmp22;
  const Scalar _tmp53 = _tmp28 + _tmp52;
  const Scalar _tmp54 = _TransformationMatrix[5] + _tmp53;
  const Scalar _tmp55 = _tmp54 - p_b(1, 0);
  const Scalar _tmp56 = _TransformationMatrix[4] + _tmp51;
  const Scalar _tmp57 = _tmp56 - p_b(0, 0);
  const Scalar _tmp58 =
      std::sqrt(Scalar(std::pow(_tmp55, Scalar(2)) + std::pow(_tmp57, Scalar(2))));
  const Scalar _tmp59 = Scalar(1.0) / (_tmp58);
  const Scalar _tmp60 = Scalar(1.0) / (_tmp57);
  const Scalar _tmp61 = _tmp58 * _tmp60;
  const Scalar _tmp62 = _tmp61 * (_tmp51 * _tmp55 * _tmp59 - _tmp53 * _tmp57 * _tmp59);
  const Scalar _tmp63 = _tmp47 * _tmp62;
  const Scalar _tmp64 = -_tmp17 * _tmp48 + _tmp19 * _tmp49 + _tmp19 * _tmp63;
  const Scalar _tmp65 = _tmp55 * _tmp60;
  const Scalar _tmp66 = -_tmp45 + _tmp46 * _tmp65;
  const Scalar _tmp67 = _tmp47 * _tmp65;
  const Scalar _tmp68 = _tmp19 * _tmp67 - _tmp48;
  const Scalar _tmp69 = Scalar(1.0) / (_tmp68);
  const Scalar _tmp70 = _tmp66 * _tmp69;
  const Scalar _tmp71 = -_tmp37 * _tmp45 + _tmp39 * _tmp46 + _tmp46 * _tmp62 - _tmp64 * _tmp70;
  const Scalar _tmp72 = Scalar(1.0) / (_tmp71);
  const Scalar _tmp73 = Scalar(1.0) * _tmp72;
  const Scalar _tmp74 = _tmp70 * _tmp73;
  const Scalar _tmp75 = _tmp19 * _tmp31 * _tmp34;
  const Scalar _tmp76 = -_tmp17 * _tmp75 + _tmp29 * _tmp35 + _tmp35 * _tmp62 - _tmp49 - _tmp63;
  const Scalar _tmp77 = std::pow(_tmp68, Scalar(-2));
  const Scalar _tmp78 = _tmp66 * _tmp77;
  const Scalar _tmp79 = _tmp35 * _tmp65 - _tmp67 - _tmp75;
  const Scalar _tmp80 = _tmp64 * _tmp79;
  const Scalar _tmp81 = -_tmp70 * _tmp76 + _tmp78 * _tmp80;
  const Scalar _tmp82 = std::pow(_tmp71, Scalar(-2));
  const Scalar _tmp83 = _tmp81 * _tmp82;
  const Scalar _tmp84 = _tmp19 * _tmp47;
  const Scalar _tmp85 = _tmp70 * _tmp84;
  const Scalar _tmp86 = Scalar(1.0) * _tmp85;
  const Scalar _tmp87 = _tmp78 * _tmp79;
  const Scalar _tmp88 = _tmp73 * _tmp84;
  const Scalar _tmp89 = Scalar(1.0) * _tmp46;
  const Scalar _tmp90 = _tmp36 + _tmp50;
  const Scalar _tmp91 = _TransformationMatrix[4] + _tmp90 - p_a(0, 0);
  const Scalar _tmp92 = _tmp38 + _tmp52;
  const Scalar _tmp93 = _TransformationMatrix[5] + _tmp92 - p_a(1, 0);
  const Scalar _tmp94 = std::pow(Scalar(std::pow(_tmp91, Scalar(2)) + std::pow(_tmp93, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp95 = _tmp93 * _tmp94;
  const Scalar _tmp96 = _tmp91 * _tmp94;
  const Scalar _tmp97 = fh1 * (_tmp90 * _tmp95 - _tmp92 * _tmp96);
  const Scalar _tmp98 = _tmp61 * _tmp97;
  const Scalar _tmp99 = _tmp65 * _tmp69;
  const Scalar _tmp100 =
      Scalar(0.20999999999999999) * _tmp24 + Scalar(0.20999999999999999) * _tmp25;
  const Scalar _tmp101 = -_tmp100;
  const Scalar _tmp102 = -Scalar(0.010999999999999999) * _tmp21 -
                         Scalar(0.010999999999999999) * _tmp9 + Scalar(-0.010999999999999999);
  const Scalar _tmp103 =
      Scalar(0.20999999999999999) * _tmp12 - Scalar(0.20999999999999999) * _tmp13;
  const Scalar _tmp104 = _tmp102 + _tmp103;
  const Scalar _tmp105 = _tmp101 + _tmp104;
  const Scalar _tmp106 = _tmp105 * _tmp47;
  const Scalar _tmp107 = _tmp106 * _tmp19;
  const Scalar _tmp108 = _tmp100 + _tmp104;
  const Scalar _tmp109 = -_tmp107 * _tmp65 + _tmp108 * _tmp48;
  const Scalar _tmp110 = _tmp105 * _tmp65;
  const Scalar _tmp111 = _tmp109 * _tmp99 + _tmp110;
  const Scalar _tmp112 = Scalar(1.0) * _tmp53;
  const Scalar _tmp113 = -_tmp112;
  const Scalar _tmp114 = Scalar(1.0) / (_tmp113 + _tmp29);
  const Scalar _tmp115 = Scalar(1.0) * _tmp51;
  const Scalar _tmp116 = _tmp115 - _tmp17;
  const Scalar _tmp117 = _tmp114 * _tmp116;
  const Scalar _tmp118 = _tmp108 * _tmp47;
  const Scalar _tmp119 = _tmp107 - _tmp118 * _tmp19;
  const Scalar _tmp120 = -_tmp105 - _tmp111 * _tmp117 + _tmp119 * _tmp99;
  const Scalar _tmp121 = _tmp102 - _tmp103;
  const Scalar _tmp122 = _tmp100 + _tmp121;
  const Scalar _tmp123 = -_tmp109 * _tmp70 - _tmp110 * _tmp46 + _tmp122 * _tmp45;
  const Scalar _tmp124 = _tmp105 * _tmp46 - _tmp117 * _tmp123 - _tmp119 * _tmp70 - _tmp122 * _tmp46;
  const Scalar _tmp125 = Scalar(1.0) / (_tmp124);
  const Scalar _tmp126 = _tmp125 * _tmp71;
  const Scalar _tmp127 = -_tmp120 * _tmp126 - _tmp62 + _tmp64 * _tmp99;
  const Scalar _tmp128 = _tmp124 * _tmp72;
  const Scalar _tmp129 = _tmp127 * _tmp128;
  const Scalar _tmp130 = _tmp120 + _tmp129;
  const Scalar _tmp131 = _tmp125 * _tmp66;
  const Scalar _tmp132 = -_tmp130 * _tmp131 - _tmp65;
  const Scalar _tmp133 = _tmp77 * _tmp79;
  const Scalar _tmp134 = _tmp133 * _tmp84;
  const Scalar _tmp135 = _tmp35 * _tmp69;
  const Scalar _tmp136 = _tmp109 * _tmp79;
  const Scalar _tmp137 = _tmp108 * _tmp75;
  const Scalar _tmp138 = _tmp106 * _tmp65 - _tmp110 * _tmp35 + _tmp137;
  const Scalar _tmp139 = _tmp136 * _tmp78 - _tmp138 * _tmp70;
  const Scalar _tmp140 = _tmp105 * _tmp35 - _tmp106 - _tmp108 * _tmp35 + _tmp118;
  const Scalar _tmp141 = _tmp119 * _tmp79;
  const Scalar _tmp142 = -_tmp117 * _tmp139 - _tmp140 * _tmp70 + _tmp141 * _tmp78;
  const Scalar _tmp143 = std::pow(_tmp124, Scalar(-2));
  const Scalar _tmp144 = _tmp142 * _tmp143;
  const Scalar _tmp145 = _tmp130 * _tmp144;
  const Scalar _tmp146 = _tmp65 * _tmp77;
  const Scalar _tmp147 = -_tmp136 * _tmp146 + _tmp138 * _tmp99;
  const Scalar _tmp148 = -_tmp117 * _tmp147 + _tmp140 * _tmp99 - _tmp141 * _tmp146;
  const Scalar _tmp149 = _tmp120 * _tmp125;
  const Scalar _tmp150 = _tmp120 * _tmp71;
  const Scalar _tmp151 = _tmp128 * (-_tmp126 * _tmp148 + _tmp144 * _tmp150 - _tmp146 * _tmp80 -
                                    _tmp149 * _tmp81 + _tmp76 * _tmp99);
  const Scalar _tmp152 = _tmp124 * _tmp83;
  const Scalar _tmp153 = _tmp127 * _tmp152;
  const Scalar _tmp154 = _tmp142 * _tmp72;
  const Scalar _tmp155 = _tmp127 * _tmp154;
  const Scalar _tmp156 = _tmp148 + _tmp151 - _tmp153 + _tmp155;
  const Scalar _tmp157 = -_tmp131 * _tmp156 + _tmp145 * _tmp66;
  const Scalar _tmp158 = _tmp47 * _tmp69;
  const Scalar _tmp159 = _tmp158 * _tmp19;
  const Scalar _tmp160 = _tmp132 * _tmp158;
  const Scalar _tmp161 = _tmp125 * _tmp46;
  const Scalar _tmp162 = _tmp96 * fh1;
  const Scalar _tmp163 = _tmp162 * _tmp61;
  const Scalar _tmp164 = _tmp112 * _tmp117 + _tmp115;
  const Scalar _tmp165 = 0;
  const Scalar _tmp166 = _tmp131 * _tmp165;
  const Scalar _tmp167 = _tmp165 * _tmp46;
  const Scalar _tmp168 = _tmp165 * _tmp85;
  const Scalar _tmp169 = Scalar(333.54000000000002) - fv1;
  const Scalar _tmp170 = _tmp169 * _tmp61;
  const Scalar _tmp171 = Scalar(1.0) * _tmp69;
  const Scalar _tmp172 = Scalar(1.0) * _tmp114;
  const Scalar _tmp173 = _tmp116 * _tmp172;
  const Scalar _tmp174 = _tmp173 * _tmp69;
  const Scalar _tmp175 = _tmp109 * _tmp174 - _tmp119 * _tmp171;
  const Scalar _tmp176 = -_tmp126 * _tmp175 - _tmp171 * _tmp64;
  const Scalar _tmp177 = _tmp128 * _tmp176;
  const Scalar _tmp178 = _tmp175 + _tmp177;
  const Scalar _tmp179 = -_tmp131 * _tmp178 + Scalar(1.0);
  const Scalar _tmp180 = _tmp158 * _tmp179;
  const Scalar _tmp181 = _tmp144 * _tmp178;
  const Scalar _tmp182 = _tmp152 * _tmp176;
  const Scalar _tmp183 = _tmp125 * _tmp175;
  const Scalar _tmp184 = Scalar(1.0) * _tmp77;
  const Scalar _tmp185 =
      -_tmp136 * _tmp173 * _tmp77 + _tmp138 * _tmp174 - _tmp140 * _tmp171 + _tmp141 * _tmp184;
  const Scalar _tmp186 = _tmp175 * _tmp71;
  const Scalar _tmp187 = _tmp128 * (-_tmp126 * _tmp185 + _tmp144 * _tmp186 - _tmp171 * _tmp76 -
                                    _tmp183 * _tmp81 + _tmp184 * _tmp80);
  const Scalar _tmp188 = _tmp154 * _tmp176;
  const Scalar _tmp189 = -_tmp182 + _tmp185 + _tmp187 + _tmp188;
  const Scalar _tmp190 = -_tmp131 * _tmp189 + _tmp181 * _tmp66;
  const Scalar _tmp191 = _tmp95 * fh1;
  const Scalar _tmp192 = _tmp191 * _tmp61;
  const Scalar _tmp193 = -_tmp163 * (-_tmp132 * _tmp134 + _tmp132 * _tmp135 - _tmp145 * _tmp46 +
                                     _tmp156 * _tmp161 + _tmp157 * _tmp159 - _tmp160) -
                         _tmp170 * (_tmp134 * _tmp166 - _tmp135 * _tmp166 - _tmp144 * _tmp167 +
                                    _tmp144 * _tmp168 + _tmp158 * _tmp166) -
                         _tmp192 * (-_tmp134 * _tmp179 + _tmp135 * _tmp179 + _tmp159 * _tmp190 +
                                    _tmp161 * _tmp189 - _tmp180 - _tmp181 * _tmp46) -
                         _tmp98 * (-_tmp35 * _tmp74 + _tmp47 * _tmp74 + _tmp83 * _tmp86 -
                                   _tmp83 * _tmp89 + _tmp87 * _tmp88);
  const Scalar _tmp194 = fh1 * (_tmp101 + _tmp121);
  const Scalar _tmp195 = -_tmp194 * _tmp95 - Scalar(40.024799999999999) * _tmp26 - _tmp92 * fv1;
  const Scalar _tmp196 = _tmp113 + _tmp39;
  const Scalar _tmp197 = _tmp117 * _tmp196;
  const Scalar _tmp198 = Scalar(1.0) / (_tmp115 - _tmp197 - _tmp37);
  const Scalar _tmp199 = Scalar(1.0) * _tmp198;
  const Scalar _tmp200 = _tmp117 * _tmp199;
  const Scalar _tmp201 = _tmp197 * _tmp199 + Scalar(1.0);
  const Scalar _tmp202 = _tmp128 * _tmp199;
  const Scalar _tmp203 = -_tmp123 * _tmp73 + _tmp196 * _tmp202;
  const Scalar _tmp204 = Scalar(1.0) * _tmp97;
  const Scalar _tmp205 = _tmp123 * _tmp125;
  const Scalar _tmp206 = _tmp196 * _tmp198;
  const Scalar _tmp207 = _tmp111 + _tmp129 * _tmp206 - _tmp130 * _tmp205;
  const Scalar _tmp208 = Scalar(1.0) * _tmp162;
  const Scalar _tmp209 = -_tmp109 * _tmp171 + _tmp177 * _tmp206 - _tmp178 * _tmp205;
  const Scalar _tmp210 = Scalar(1.0) * _tmp191;
  const Scalar _tmp211 = _tmp114 * _tmp196;
  const Scalar _tmp212 = Scalar(40.024799999999999) * _tmp14 + _tmp194 * _tmp96 + _tmp90 * fv1;
  const Scalar _tmp213 = _tmp164 * _tmp198;
  const Scalar _tmp214 = _tmp113 - _tmp165 * _tmp205 - _tmp196 * _tmp213;
  const Scalar _tmp215 =
      Scalar(1.0) * _tmp169 * (-_tmp164 * _tmp199 - _tmp172 * _tmp214 + Scalar(1.0)) +
      Scalar(1.0) * _tmp195 * (-_tmp172 * _tmp201 + _tmp200) +
      _tmp204 * (-_tmp172 * _tmp203 + _tmp202) + _tmp208 * (_tmp129 * _tmp199 - _tmp172 * _tmp207) +
      _tmp210 * (-_tmp172 * _tmp209 + _tmp177 * _tmp199) +
      Scalar(1.0) * _tmp212 * (_tmp199 * _tmp211 - _tmp199);
  const Scalar _tmp216 = _tmp125 * _tmp165;
  const Scalar _tmp217 = _tmp125 * _tmp178;
  const Scalar _tmp218 = _tmp125 * _tmp130;
  const Scalar _tmp219 = -_tmp163 * (_tmp160 * _tmp19 + _tmp218 * _tmp46 + Scalar(1.0)) -
                         _tmp170 * (-_tmp159 * _tmp166 + _tmp216 * _tmp46) -
                         _tmp192 * (_tmp180 * _tmp19 + _tmp217 * _tmp46) -
                         _tmp98 * (_tmp46 * _tmp73 - _tmp74 * _tmp84);
  const Scalar _tmp220 = std::pow(_tmp219, Scalar(-2));
  const Scalar _tmp221 = _tmp215 * _tmp220;
  const Scalar _tmp222 = Scalar(1.0) / (_tmp219);
  const Scalar _tmp223 = _tmp152 * _tmp199;
  const Scalar _tmp224 = _tmp154 * _tmp199;
  const Scalar _tmp225 = Scalar(1.0) * _tmp123;
  const Scalar _tmp226 =
      -_tmp139 * _tmp73 - _tmp196 * _tmp223 + _tmp196 * _tmp224 + _tmp225 * _tmp83;
  const Scalar _tmp227 = _tmp123 * _tmp145 - _tmp139 * _tmp218 + _tmp147 + _tmp151 * _tmp206 -
                         _tmp153 * _tmp206 + _tmp155 * _tmp206 - _tmp156 * _tmp205;
  const Scalar _tmp228 = _tmp123 * _tmp181 + _tmp136 * _tmp184 - _tmp138 * _tmp171 -
                         _tmp139 * _tmp217 - _tmp182 * _tmp206 + _tmp187 * _tmp206 +
                         _tmp188 * _tmp206 - _tmp189 * _tmp205;
  const Scalar _tmp229 = _tmp123 * _tmp165;
  const Scalar _tmp230 = -_tmp139 * _tmp216 + _tmp144 * _tmp229;
  const Scalar _tmp231 = _tmp169 * _tmp172;
  const Scalar _tmp232 =
      -_tmp193 * _tmp221 + _tmp222 * (_tmp204 * (-_tmp172 * _tmp226 - _tmp223 + _tmp224) +
                                      _tmp208 * (_tmp151 * _tmp199 - _tmp153 * _tmp199 +
                                                 _tmp155 * _tmp199 - _tmp172 * _tmp227) +
                                      _tmp210 * (-_tmp172 * _tmp228 - _tmp182 * _tmp199 +
                                                 _tmp187 * _tmp199 + _tmp188 * _tmp199) -
                                      _tmp230 * _tmp231);
  const Scalar _tmp233 =
      std::pow(Scalar(std::pow(_tmp215, Scalar(2)) * _tmp220 + 1), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp234 = Scalar(1.4083112389913199) * _tmp219;
  const Scalar _tmp235 = _tmp233 * _tmp234;
  const Scalar _tmp236 = std::asinh(_tmp215 * _tmp222);
  const Scalar _tmp237 = Scalar(1.4083112389913199) * _tmp236;
  const Scalar _tmp238 = Scalar(0.71007031138673404) * _tmp222;
  const Scalar _tmp239 =
      -_tmp234 * _tmp236 - std::sqrt(Scalar(std::pow(Scalar(-_tmp54 + p_b(1, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp56 + p_b(0, 0)), Scalar(2))));
  const Scalar _tmp240 = Scalar(0.71007031138673404) * _tmp220 * _tmp239;
  const Scalar _tmp241 = _tmp238 * _tmp239;
  const Scalar _tmp242 = std::cosh(_tmp241);
  const Scalar _tmp243 = Scalar(1.0) * _tmp236;
  const Scalar _tmp244 = Scalar(1.0) * _tmp233 * std::cosh(_tmp243);
  const Scalar _tmp245 = -Scalar(1.4083112389913199) * std::sinh(_tmp241) -
                         Scalar(1.4083112389913199) * std::sinh(_tmp243);
  const Scalar _tmp246 = -_tmp18 + p_c(0, 0);
  const Scalar _tmp247 = -_tmp30 + p_c(1, 0);
  const Scalar _tmp248 =
      std::sqrt(Scalar(std::pow(_tmp246, Scalar(2)) + std::pow(_tmp247, Scalar(2))));
  const Scalar _tmp249 = Scalar(1.0) / (_tmp248);
  const Scalar _tmp250 = _tmp73 * _tmp97;
  const Scalar _tmp251 = _tmp204 * _tmp83;
  const Scalar _tmp252 = _tmp191 * _tmp69;
  const Scalar _tmp253 = _tmp165 * _tmp169;
  const Scalar _tmp254 = _tmp144 * _tmp253;
  const Scalar _tmp255 = _tmp131 * _tmp253;
  const Scalar _tmp256 = _tmp162 * _tmp69;
  const Scalar _tmp257 = -_tmp132 * _tmp133 * _tmp162 - _tmp133 * _tmp179 * _tmp191 +
                         _tmp133 * _tmp255 + _tmp157 * _tmp256 + _tmp190 * _tmp252 +
                         _tmp250 * _tmp87 + _tmp251 * _tmp70 + _tmp254 * _tmp70;
  const Scalar _tmp258 = _tmp114 * _tmp97;
  const Scalar _tmp259 = _tmp114 * _tmp169;
  const Scalar _tmp260 = _tmp114 * _tmp162;
  const Scalar _tmp261 = _tmp199 * _tmp212;
  const Scalar _tmp262 = _tmp114 * _tmp191;
  const Scalar _tmp263 = _tmp114 * _tmp195 * _tmp201 + _tmp203 * _tmp258 + _tmp207 * _tmp260 +
                         _tmp209 * _tmp262 - _tmp211 * _tmp261 + _tmp214 * _tmp259;
  const Scalar _tmp264 =
      _tmp132 * _tmp256 + _tmp179 * _tmp252 - _tmp250 * _tmp70 - _tmp255 * _tmp69;
  const Scalar _tmp265 = Scalar(1.0) / (_tmp264);
  const Scalar _tmp266 = std::asinh(_tmp263 * _tmp265);
  const Scalar _tmp267 = Scalar(1.4083112389913199) * _tmp266;
  const Scalar _tmp268 = std::pow(_tmp264, Scalar(-2));
  const Scalar _tmp269 = _tmp263 * _tmp268;
  const Scalar _tmp270 = -_tmp257 * _tmp269 + _tmp265 * (_tmp226 * _tmp258 + _tmp227 * _tmp260 +
                                                         _tmp228 * _tmp262 + _tmp230 * _tmp259);
  const Scalar _tmp271 =
      std::pow(Scalar(std::pow(_tmp263, Scalar(2)) * _tmp268 + 1), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp272 = Scalar(1.4083112389913199) * _tmp264;
  const Scalar _tmp273 = _tmp271 * _tmp272;
  const Scalar _tmp274 = Scalar(0.71007031138673404) * _tmp265;
  const Scalar _tmp275 = -_tmp248 - _tmp264 * _tmp267;
  const Scalar _tmp276 = Scalar(0.71007031138673404) * _tmp268 * _tmp275;
  const Scalar _tmp277 = _tmp274 * _tmp275;
  const Scalar _tmp278 = std::cosh(_tmp277);
  const Scalar _tmp279 = Scalar(1.0) * _tmp266;
  const Scalar _tmp280 = Scalar(1.0) * _tmp271 * std::cosh(_tmp279);
  const Scalar _tmp281 = -Scalar(1.4083112389913199) * std::sinh(_tmp277) -
                         Scalar(1.4083112389913199) * std::sinh(_tmp279);
  const Scalar _tmp282 = _tmp162 * _tmp218 + _tmp169 * _tmp216 + _tmp191 * _tmp217 + _tmp250;
  const Scalar _tmp283 = Scalar(1.4083112389913199) * _tmp282;
  const Scalar _tmp284 = _tmp191 * _tmp198;
  const Scalar _tmp285 = _tmp162 * _tmp198;
  const Scalar _tmp286 = -_tmp129 * _tmp285 + _tmp169 * _tmp213 - _tmp177 * _tmp284 -
                         _tmp195 * _tmp200 - _tmp202 * _tmp97 + _tmp261;
  const Scalar _tmp287 = std::pow(_tmp282, Scalar(-2));
  const Scalar _tmp288 =
      std::pow(Scalar(std::pow(_tmp286, Scalar(2)) * _tmp287 + 1), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp289 = _tmp125 * _tmp162;
  const Scalar _tmp290 = _tmp125 * _tmp191;
  const Scalar _tmp291 = -_tmp145 * _tmp162 + _tmp156 * _tmp289 - _tmp181 * _tmp191 +
                         _tmp189 * _tmp290 - _tmp251 - _tmp254;
  const Scalar _tmp292 = _tmp286 * _tmp287;
  const Scalar _tmp293 = Scalar(1.0) / (_tmp282);
  const Scalar _tmp294 =
      _tmp288 *
      (-_tmp291 * _tmp292 +
       _tmp293 * (-_tmp151 * _tmp285 + _tmp153 * _tmp285 - _tmp155 * _tmp285 + _tmp182 * _tmp284 -
                  _tmp187 * _tmp284 - _tmp188 * _tmp284 + _tmp223 * _tmp97 - _tmp224 * _tmp97));
  const Scalar _tmp295 = std::asinh(_tmp286 * _tmp293);
  const Scalar _tmp296 = Scalar(1.4083112389913199) * _tmp295;
  const Scalar _tmp297 = Scalar(0.71007031138673404) * _tmp293;
  const Scalar _tmp298 =
      -_tmp283 * _tmp295 - std::sqrt(Scalar(std::pow(Scalar(-_tmp40 + p_d(1, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp42 + p_d(0, 0)), Scalar(2))));
  const Scalar _tmp299 = Scalar(0.71007031138673404) * _tmp287 * _tmp298;
  const Scalar _tmp300 = _tmp297 * _tmp298;
  const Scalar _tmp301 = std::cosh(_tmp300);
  const Scalar _tmp302 = Scalar(1.0) * _tmp295;
  const Scalar _tmp303 = Scalar(1.0) * std::cosh(_tmp302);
  const Scalar _tmp304 = -Scalar(1.4083112389913199) * std::sinh(_tmp300) -
                         Scalar(1.4083112389913199) * std::sinh(_tmp302);
  const Scalar _tmp305 = _tmp32 * _tmp34;
  const Scalar _tmp306 = _tmp108 * _tmp305 - _tmp110 * _tmp75 - _tmp118;
  const Scalar _tmp307 = -_tmp305 + _tmp47 + _tmp65 * _tmp75;
  const Scalar _tmp308 = _tmp307 * _tmp78;
  const Scalar _tmp309 = _tmp109 * _tmp308 - _tmp306 * _tmp70;
  const Scalar _tmp310 = _tmp105 * _tmp75 - _tmp137;
  const Scalar _tmp311 = -_tmp117 * _tmp309 + _tmp119 * _tmp308 - _tmp310 * _tmp70;
  const Scalar _tmp312 = _tmp143 * _tmp311;
  const Scalar _tmp313 = _tmp178 * _tmp312;
  const Scalar _tmp314 = _tmp69 * _tmp75;
  const Scalar _tmp315 = -_tmp17 * _tmp305 + _tmp17 * _tmp47 + _tmp29 * _tmp75 + _tmp62 * _tmp75;
  const Scalar _tmp316 = _tmp308 * _tmp64 - _tmp315 * _tmp70;
  const Scalar _tmp317 = _tmp316 * _tmp82;
  const Scalar _tmp318 = _tmp124 * _tmp317;
  const Scalar _tmp319 = _tmp176 * _tmp318;
  const Scalar _tmp320 = _tmp307 * _tmp77;
  const Scalar _tmp321 = _tmp184 * _tmp307;
  const Scalar _tmp322 =
      -_tmp109 * _tmp173 * _tmp320 + _tmp119 * _tmp321 - _tmp171 * _tmp310 + _tmp174 * _tmp306;
  const Scalar _tmp323 = _tmp128 * (-_tmp126 * _tmp322 - _tmp171 * _tmp315 - _tmp183 * _tmp316 +
                                    _tmp186 * _tmp312 + _tmp321 * _tmp64);
  const Scalar _tmp324 = _tmp311 * _tmp72;
  const Scalar _tmp325 = _tmp176 * _tmp324;
  const Scalar _tmp326 = -_tmp319 + _tmp322 + _tmp323 + _tmp325;
  const Scalar _tmp327 = -_tmp131 * _tmp326 + _tmp313 * _tmp66;
  const Scalar _tmp328 = _tmp179 * _tmp320;
  const Scalar _tmp329 = _tmp132 * _tmp320;
  const Scalar _tmp330 = _tmp127 * _tmp318;
  const Scalar _tmp331 = _tmp320 * _tmp65;
  const Scalar _tmp332 = -_tmp109 * _tmp331 + _tmp306 * _tmp99;
  const Scalar _tmp333 = -_tmp117 * _tmp332 - _tmp119 * _tmp331 + _tmp310 * _tmp99;
  const Scalar _tmp334 = _tmp128 * (-_tmp126 * _tmp333 - _tmp149 * _tmp316 + _tmp150 * _tmp312 +
                                    _tmp315 * _tmp99 - _tmp331 * _tmp64);
  const Scalar _tmp335 = _tmp127 * _tmp324;
  const Scalar _tmp336 = -_tmp330 + _tmp333 + _tmp334 + _tmp335;
  const Scalar _tmp337 = _tmp130 * _tmp312;
  const Scalar _tmp338 = -_tmp131 * _tmp336 + _tmp337 * _tmp66;
  const Scalar _tmp339 =
      -_tmp163 * (_tmp132 * _tmp314 + _tmp159 * _tmp338 + _tmp161 * _tmp336 - _tmp329 * _tmp84 -
                  _tmp337 * _tmp46) -
      _tmp170 * (-_tmp166 * _tmp314 + _tmp166 * _tmp320 * _tmp84 - _tmp167 * _tmp312 +
                 _tmp168 * _tmp312) -
      _tmp192 * (_tmp159 * _tmp327 + _tmp161 * _tmp326 + _tmp179 * _tmp314 - _tmp313 * _tmp46 -
                 _tmp328 * _tmp84) -
      _tmp98 * (_tmp308 * _tmp88 + _tmp317 * _tmp86 - _tmp317 * _tmp89 - _tmp74 * _tmp75);
  const Scalar _tmp340 = _tmp199 * _tmp318;
  const Scalar _tmp341 = _tmp199 * _tmp324;
  const Scalar _tmp342 =
      -_tmp196 * _tmp340 + _tmp196 * _tmp341 + _tmp225 * _tmp317 - _tmp309 * _tmp73;
  const Scalar _tmp343 = -_tmp216 * _tmp309 + _tmp229 * _tmp312;
  const Scalar _tmp344 = _tmp123 * _tmp337 - _tmp205 * _tmp336 - _tmp206 * _tmp330 +
                         _tmp206 * _tmp334 + _tmp206 * _tmp335 - _tmp218 * _tmp309 + _tmp332;
  const Scalar _tmp345 = _tmp109 * _tmp321 + _tmp123 * _tmp313 - _tmp171 * _tmp306 -
                         _tmp205 * _tmp326 - _tmp206 * _tmp319 + _tmp206 * _tmp323 +
                         _tmp206 * _tmp325 - _tmp217 * _tmp309;
  const Scalar _tmp346 =
      -_tmp221 * _tmp339 + _tmp222 * (_tmp204 * (-_tmp172 * _tmp342 - _tmp340 + _tmp341) +
                                      _tmp208 * (-_tmp172 * _tmp344 - _tmp199 * _tmp330 +
                                                 _tmp199 * _tmp334 + _tmp199 * _tmp335) +
                                      _tmp210 * (-_tmp172 * _tmp345 - _tmp199 * _tmp319 +
                                                 _tmp199 * _tmp323 + _tmp199 * _tmp325) -
                                      _tmp231 * _tmp343);
  const Scalar _tmp347 = _tmp204 * _tmp317;
  const Scalar _tmp348 = _tmp253 * _tmp312;
  const Scalar _tmp349 = -_tmp162 * _tmp329 - _tmp191 * _tmp328 + _tmp250 * _tmp308 +
                         _tmp252 * _tmp327 + _tmp255 * _tmp320 + _tmp256 * _tmp338 +
                         _tmp347 * _tmp70 + _tmp348 * _tmp70;
  const Scalar _tmp350 =
      _tmp265 * (_tmp258 * _tmp342 + _tmp259 * _tmp343 + _tmp260 * _tmp344 + _tmp262 * _tmp345) -
      _tmp269 * _tmp349;
  const Scalar _tmp351 = -_tmp162 * _tmp337 - _tmp191 * _tmp313 + _tmp289 * _tmp336 +
                         _tmp290 * _tmp326 - _tmp347 - _tmp348;
  const Scalar _tmp352 =
      _tmp288 *
      (-_tmp292 * _tmp351 +
       _tmp293 * (_tmp284 * _tmp319 - _tmp284 * _tmp323 - _tmp284 * _tmp325 + _tmp285 * _tmp330 -
                  _tmp285 * _tmp334 - _tmp285 * _tmp335 + _tmp340 * _tmp97 - _tmp341 * _tmp97));

  // Output terms (1)
  Eigen::Matrix<Scalar, 4, 3> _res;

  _res(0, 0) = 0;
  _res(1, 0) =
      _tmp193 * _tmp245 +
      _tmp234 *
          (-_tmp232 * _tmp244 -
           _tmp242 * (-_tmp193 * _tmp240 + _tmp238 * (-_tmp193 * _tmp237 - _tmp232 * _tmp235)));
  _res(2, 0) =
      _tmp257 * _tmp281 +
      _tmp272 * (-_tmp270 * _tmp280 -
                 _tmp278 * (-_tmp257 * _tmp276 + _tmp274 * (-_tmp246 * _tmp249 - _tmp257 * _tmp267 -
                                                            _tmp270 * _tmp273)));
  _res(3, 0) =
      _tmp283 *
          (-_tmp294 * _tmp303 -
           _tmp301 * (-_tmp291 * _tmp299 + _tmp297 * (-_tmp283 * _tmp294 - _tmp291 * _tmp296))) +
      _tmp291 * _tmp304;
  _res(0, 1) = 0;
  _res(1, 1) =
      _tmp234 *
          (-_tmp242 * (_tmp238 * (-_tmp235 * _tmp346 - _tmp237 * _tmp339) - _tmp240 * _tmp339) -
           _tmp244 * _tmp346) +
      _tmp245 * _tmp339;
  _res(2, 1) =
      _tmp272 *
          (-_tmp278 * (_tmp274 * (-_tmp247 * _tmp249 - _tmp267 * _tmp349 - _tmp273 * _tmp350) -
                       _tmp276 * _tmp349) -
           _tmp280 * _tmp350) +
      _tmp281 * _tmp349;
  _res(3, 1) =
      _tmp283 *
          (-_tmp301 * (_tmp297 * (-_tmp283 * _tmp352 - _tmp296 * _tmp351) - _tmp299 * _tmp351) -
           _tmp303 * _tmp352) +
      _tmp304 * _tmp351;
  _res(0, 2) = 0;
  _res(1, 2) = 0;
  _res(2, 2) = 0;
  _res(3, 2) = 0;

  return _res;
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym

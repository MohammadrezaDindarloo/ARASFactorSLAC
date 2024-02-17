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
  // Total ops: 1138

  // Unused inputs
  (void)encoder;
  (void)offset;
  (void)epsilon;

  // Input arrays
  const Eigen::Matrix<Scalar, 4, 1>& _DeltaRot = DeltaRot.Data();
  const Eigen::Matrix<Scalar, 7, 1>& _TransformationMatrix = TransformationMatrix.Data();

  // Intermediate terms (358)
  const Scalar _tmp0 =
      _DeltaRot[0] * _TransformationMatrix[2] + _DeltaRot[1] * _TransformationMatrix[3] -
      _DeltaRot[2] * _TransformationMatrix[0] + _DeltaRot[3] * _TransformationMatrix[1];
  const Scalar _tmp1 = -2 * std::pow(_tmp0, Scalar(2));
  const Scalar _tmp2 =
      -_DeltaRot[0] * _TransformationMatrix[1] + _DeltaRot[1] * _TransformationMatrix[0] +
      _DeltaRot[2] * _TransformationMatrix[3] + _DeltaRot[3] * _TransformationMatrix[2];
  const Scalar _tmp3 = 1 - 2 * std::pow(_tmp2, Scalar(2));
  const Scalar _tmp4 = Scalar(0.20999999999999999) * _tmp1 + Scalar(0.20999999999999999) * _tmp3;
  const Scalar _tmp5 =
      _DeltaRot[0] * _TransformationMatrix[3] - _DeltaRot[1] * _TransformationMatrix[2] +
      _DeltaRot[2] * _TransformationMatrix[1] + _DeltaRot[3] * _TransformationMatrix[0];
  const Scalar _tmp6 = 2 * _tmp5;
  const Scalar _tmp7 = _tmp2 * _tmp6;
  const Scalar _tmp8 =
      -2 * _DeltaRot[0] * _TransformationMatrix[0] - 2 * _DeltaRot[1] * _TransformationMatrix[1] -
      2 * _DeltaRot[2] * _TransformationMatrix[2] + 2 * _DeltaRot[3] * _TransformationMatrix[3];
  const Scalar _tmp9 = _tmp0 * _tmp8;
  const Scalar _tmp10 = _tmp7 + _tmp9;
  const Scalar _tmp11 = -Scalar(0.010999999999999999) * _tmp10;
  const Scalar _tmp12 = _tmp0 * _tmp6;
  const Scalar _tmp13 = _tmp2 * _tmp8;
  const Scalar _tmp14 = Scalar(0.20999999999999999) * _tmp12 - Scalar(0.20999999999999999) * _tmp13;
  const Scalar _tmp15 = _tmp11 - _tmp14;
  const Scalar _tmp16 = _tmp15 + _tmp4;
  const Scalar _tmp17 = _TransformationMatrix[4] + _tmp16;
  const Scalar _tmp18 = Scalar(0.20999999999999999) * _tmp12 + Scalar(0.20999999999999999) * _tmp13;
  const Scalar _tmp19 = 2 * _tmp0 * _tmp2;
  const Scalar _tmp20 = _tmp5 * _tmp8;
  const Scalar _tmp21 = _tmp19 - _tmp20;
  const Scalar _tmp22 = -Scalar(0.010999999999999999) * _tmp21;
  const Scalar _tmp23 = -2 * std::pow(_tmp5, Scalar(2));
  const Scalar _tmp24 = Scalar(0.20999999999999999) * _tmp23 + Scalar(0.20999999999999999) * _tmp3;
  const Scalar _tmp25 = _tmp22 - _tmp24;
  const Scalar _tmp26 = _tmp18 + _tmp25;
  const Scalar _tmp27 = _TransformationMatrix[5] + _tmp26;
  const Scalar _tmp28 = _tmp11 + _tmp14;
  const Scalar _tmp29 = _tmp28 + _tmp4;
  const Scalar _tmp30 = _TransformationMatrix[4] + _tmp29;
  const Scalar _tmp31 = _tmp30 - p_c(0, 0);
  const Scalar _tmp32 = std::pow(_tmp31, Scalar(2));
  const Scalar _tmp33 = _tmp22 + _tmp24;
  const Scalar _tmp34 = _tmp18 + _tmp33;
  const Scalar _tmp35 = _TransformationMatrix[5] + _tmp34;
  const Scalar _tmp36 = _tmp35 - p_c(1, 0);
  const Scalar _tmp37 = std::pow(_tmp36, Scalar(2));
  const Scalar _tmp38 = _tmp32 + _tmp37;
  const Scalar _tmp39 = std::pow(_tmp38, Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp40 = _tmp31 * _tmp39;
  const Scalar _tmp41 = -_tmp18;
  const Scalar _tmp42 = _tmp33 + _tmp41;
  const Scalar _tmp43 = _TransformationMatrix[5] + _tmp42;
  const Scalar _tmp44 = _tmp43 - p_d(1, 0);
  const Scalar _tmp45 = -_tmp4;
  const Scalar _tmp46 = _tmp28 + _tmp45;
  const Scalar _tmp47 = _TransformationMatrix[4] + _tmp46;
  const Scalar _tmp48 = _tmp47 - p_d(0, 0);
  const Scalar _tmp49 = std::pow(Scalar(std::pow(_tmp44, Scalar(2)) + std::pow(_tmp48, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp50 = _tmp44 * _tmp49;
  const Scalar _tmp51 = _tmp27 - p_b(1, 0);
  const Scalar _tmp52 = _tmp17 - p_b(0, 0);
  const Scalar _tmp53 = Scalar(1.0) / (_tmp52);
  const Scalar _tmp54 = _tmp51 * _tmp53;
  const Scalar _tmp55 = _tmp48 * _tmp49;
  const Scalar _tmp56 = -_tmp50 + _tmp54 * _tmp55;
  const Scalar _tmp57 = _tmp39 * _tmp54;
  const Scalar _tmp58 = _tmp31 * _tmp57 - _tmp36 * _tmp39;
  const Scalar _tmp59 = Scalar(1.0) / (_tmp58);
  const Scalar _tmp60 = _tmp56 * _tmp59;
  const Scalar _tmp61 = Scalar(0.20999999999999999) * _tmp19 + Scalar(0.20999999999999999) * _tmp20;
  const Scalar _tmp62 = -Scalar(0.010999999999999999) * _tmp1 -
                        Scalar(0.010999999999999999) * _tmp23 + Scalar(-0.010999999999999999);
  const Scalar _tmp63 = Scalar(0.20999999999999999) * _tmp7 - Scalar(0.20999999999999999) * _tmp9;
  const Scalar _tmp64 = _tmp62 - _tmp63;
  const Scalar _tmp65 = _tmp61 + _tmp64;
  const Scalar _tmp66 = -_tmp61;
  const Scalar _tmp67 = _tmp62 + _tmp63;
  const Scalar _tmp68 = _tmp66 + _tmp67;
  const Scalar _tmp69 = _tmp61 + _tmp67;
  const Scalar _tmp70 = _tmp39 * _tmp69;
  const Scalar _tmp71 = _tmp39 * _tmp68;
  const Scalar _tmp72 = _tmp31 * _tmp71;
  const Scalar _tmp73 = -_tmp31 * _tmp70 + _tmp72;
  const Scalar _tmp74 = _tmp36 * _tmp70 - _tmp54 * _tmp72;
  const Scalar _tmp75 = _tmp54 * _tmp68;
  const Scalar _tmp76 = _tmp50 * _tmp65 - _tmp55 * _tmp75 - _tmp60 * _tmp74;
  const Scalar _tmp77 = Scalar(1.0) * _tmp26;
  const Scalar _tmp78 = -_tmp77;
  const Scalar _tmp79 = Scalar(1.0) / (_tmp34 + _tmp78);
  const Scalar _tmp80 = Scalar(1.0) * _tmp16;
  const Scalar _tmp81 = -_tmp29 + _tmp80;
  const Scalar _tmp82 = _tmp79 * _tmp81;
  const Scalar _tmp83 = -_tmp55 * _tmp65 + _tmp55 * _tmp68 - _tmp60 * _tmp73 - _tmp76 * _tmp82;
  const Scalar _tmp84 = Scalar(1.0) / (_tmp83);
  const Scalar _tmp85 = _tmp77 * _tmp82 + _tmp80;
  const Scalar _tmp86 = 0;
  const Scalar _tmp87 = _tmp84 * _tmp86;
  const Scalar _tmp88 = _tmp60 * _tmp87;
  const Scalar _tmp89 = _tmp55 * _tmp84;
  const Scalar _tmp90 = Scalar(43.164000000000001) - fv1;
  const Scalar _tmp91 =
      std::sqrt(Scalar(std::pow(_tmp51, Scalar(2)) + std::pow(_tmp52, Scalar(2))));
  const Scalar _tmp92 = _tmp53 * _tmp91;
  const Scalar _tmp93 = _tmp90 * _tmp92;
  const Scalar _tmp94 = Scalar(1.0) * _tmp59;
  const Scalar _tmp95 = Scalar(1.0) * _tmp79;
  const Scalar _tmp96 = _tmp81 * _tmp95;
  const Scalar _tmp97 = _tmp59 * _tmp96;
  const Scalar _tmp98 = -_tmp73 * _tmp94 + _tmp74 * _tmp97;
  const Scalar _tmp99 = _tmp29 * _tmp39;
  const Scalar _tmp100 = _tmp34 * _tmp39;
  const Scalar _tmp101 = Scalar(1.0) / (_tmp91);
  const Scalar _tmp102 = _tmp92 * (_tmp101 * _tmp16 * _tmp51 - _tmp101 * _tmp26 * _tmp52);
  const Scalar _tmp103 = _tmp102 * _tmp39;
  const Scalar _tmp104 = _tmp100 * _tmp31 + _tmp103 * _tmp31 - _tmp36 * _tmp99;
  const Scalar _tmp105 = _tmp102 * _tmp55 - _tmp104 * _tmp60 + _tmp42 * _tmp55 - _tmp46 * _tmp50;
  const Scalar _tmp106 = _tmp105 * _tmp84;
  const Scalar _tmp107 = -_tmp104 * _tmp94 - _tmp106 * _tmp98;
  const Scalar _tmp108 = Scalar(1.0) / (_tmp105);
  const Scalar _tmp109 = _tmp108 * _tmp83;
  const Scalar _tmp110 = _tmp107 * _tmp109;
  const Scalar _tmp111 = _tmp110 + _tmp98;
  const Scalar _tmp112 = _tmp56 * _tmp84;
  const Scalar _tmp113 = -_tmp111 * _tmp112 + Scalar(1.0);
  const Scalar _tmp114 = _tmp39 * _tmp59;
  const Scalar _tmp115 = _tmp113 * _tmp114;
  const Scalar _tmp116 = _tmp15 + _tmp45;
  const Scalar _tmp117 = _TransformationMatrix[4] + _tmp116 - p_a(0, 0);
  const Scalar _tmp118 = _tmp25 + _tmp41;
  const Scalar _tmp119 = _TransformationMatrix[5] + _tmp118 - p_a(1, 0);
  const Scalar _tmp120 =
      std::pow(Scalar(std::pow(_tmp117, Scalar(2)) + std::pow(_tmp119, Scalar(2))),
               Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp121 = _tmp119 * _tmp120;
  const Scalar _tmp122 = _tmp121 * fh1;
  const Scalar _tmp123 = _tmp122 * _tmp92;
  const Scalar _tmp124 = _tmp54 * _tmp59;
  const Scalar _tmp125 = _tmp124 * _tmp74 + _tmp75;
  const Scalar _tmp126 = _tmp124 * _tmp73 - _tmp125 * _tmp82 - _tmp68;
  const Scalar _tmp127 = _tmp126 * _tmp84;
  const Scalar _tmp128 = -_tmp102 + _tmp104 * _tmp124 - _tmp105 * _tmp127;
  const Scalar _tmp129 = _tmp108 * _tmp128;
  const Scalar _tmp130 = _tmp129 * _tmp83;
  const Scalar _tmp131 = _tmp126 + _tmp130;
  const Scalar _tmp132 = -_tmp112 * _tmp131 - _tmp54;
  const Scalar _tmp133 = _tmp114 * _tmp132;
  const Scalar _tmp134 = _tmp117 * _tmp120;
  const Scalar _tmp135 = _tmp134 * fh1;
  const Scalar _tmp136 = _tmp135 * _tmp92;
  const Scalar _tmp137 = Scalar(1.0) * _tmp108;
  const Scalar _tmp138 = _tmp137 * _tmp60;
  const Scalar _tmp139 = fh1 * (_tmp116 * _tmp121 - _tmp118 * _tmp134);
  const Scalar _tmp140 = _tmp139 * _tmp92;
  const Scalar _tmp141 = -_tmp123 * (_tmp111 * _tmp89 + _tmp115 * _tmp31) -
                         _tmp136 * (_tmp131 * _tmp89 + _tmp133 * _tmp31 + Scalar(1.0)) -
                         _tmp140 * (_tmp137 * _tmp55 - _tmp138 * _tmp40) -
                         _tmp93 * (-_tmp40 * _tmp88 + _tmp86 * _tmp89);
  const Scalar _tmp142 = Scalar(1.0) / (_tmp141);
  const Scalar _tmp143 = _tmp42 + _tmp78;
  const Scalar _tmp144 = _tmp143 * _tmp82;
  const Scalar _tmp145 = Scalar(1.0) / (-_tmp144 - _tmp46 + _tmp80);
  const Scalar _tmp146 = Scalar(1.0) * _tmp145;
  const Scalar _tmp147 = _tmp109 * _tmp146;
  const Scalar _tmp148 = -_tmp137 * _tmp76 + _tmp143 * _tmp147;
  const Scalar _tmp149 = Scalar(1.0) * _tmp139;
  const Scalar _tmp150 = _tmp76 * _tmp84;
  const Scalar _tmp151 = _tmp143 * _tmp145;
  const Scalar _tmp152 = _tmp125 + _tmp130 * _tmp151 - _tmp131 * _tmp150;
  const Scalar _tmp153 = Scalar(1.0) * _tmp135;
  const Scalar _tmp154 = fh1 * (_tmp64 + _tmp66);
  const Scalar _tmp155 = -_tmp118 * fv1 - _tmp121 * _tmp154 - Scalar(5.1796800000000003) * _tmp21;
  const Scalar _tmp156 = _tmp146 * _tmp82;
  const Scalar _tmp157 = _tmp144 * _tmp146 + Scalar(1.0);
  const Scalar _tmp158 = _tmp110 * _tmp151 - _tmp111 * _tmp150 - _tmp74 * _tmp94;
  const Scalar _tmp159 = Scalar(1.0) * _tmp122;
  const Scalar _tmp160 = _tmp145 * _tmp85;
  const Scalar _tmp161 = -_tmp143 * _tmp160 - _tmp150 * _tmp86 + _tmp78;
  const Scalar _tmp162 = _tmp143 * _tmp79;
  const Scalar _tmp163 = Scalar(5.1796800000000003) * _tmp10 + _tmp116 * fv1 + _tmp134 * _tmp154;
  const Scalar _tmp164 =
      _tmp149 * (_tmp147 - _tmp148 * _tmp95) + _tmp153 * (_tmp130 * _tmp146 - _tmp152 * _tmp95) +
      Scalar(1.0) * _tmp155 * (_tmp156 - _tmp157 * _tmp95) +
      _tmp159 * (_tmp110 * _tmp146 - _tmp158 * _tmp95) +
      Scalar(1.0) * _tmp163 * (_tmp146 * _tmp162 - _tmp146) +
      Scalar(1.0) * _tmp90 * (-_tmp146 * _tmp85 - _tmp161 * _tmp95 + Scalar(1.0));
  const Scalar _tmp165 = std::asinh(_tmp142 * _tmp164);
  const Scalar _tmp166 = Scalar(9.6622558468725703) * _tmp141;
  const Scalar _tmp167 =
      -_tmp165 * _tmp166 - std::sqrt(Scalar(std::pow(Scalar(-_tmp17 + p_b(0, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp27 + p_b(1, 0)), Scalar(2))));
  const Scalar _tmp168 = Scalar(0.1034955) * _tmp142;
  const Scalar _tmp169 = _tmp167 * _tmp168;
  const Scalar _tmp170 = Scalar(1.0) * _tmp165;
  const Scalar _tmp171 = -std::sinh(_tmp169) - std::sinh(_tmp170);
  const Scalar _tmp172 = std::pow(_tmp58, Scalar(-2));
  const Scalar _tmp173 = std::pow(_tmp38, Scalar(Scalar(-3) / Scalar(2)));
  const Scalar _tmp174 = _tmp173 * _tmp31 * _tmp36;
  const Scalar _tmp175 = _tmp173 * _tmp32;
  const Scalar _tmp176 = _tmp172 * (-_tmp174 + _tmp175 * _tmp54 - _tmp57);
  const Scalar _tmp177 = _tmp112 * _tmp40 * _tmp86;
  const Scalar _tmp178 = _tmp176 * _tmp56;
  const Scalar _tmp179 = _tmp174 * _tmp69;
  const Scalar _tmp180 = -_tmp175 * _tmp75 + _tmp179 + _tmp54 * _tmp71;
  const Scalar _tmp181 = _tmp178 * _tmp74 - _tmp180 * _tmp60;
  const Scalar _tmp182 = _tmp175 * _tmp68 - _tmp175 * _tmp69 + _tmp70 - _tmp71;
  const Scalar _tmp183 = _tmp178 * _tmp73 - _tmp181 * _tmp82 - _tmp182 * _tmp60;
  const Scalar _tmp184 = std::pow(_tmp83, Scalar(-2));
  const Scalar _tmp185 = _tmp183 * _tmp184;
  const Scalar _tmp186 = _tmp185 * _tmp55;
  const Scalar _tmp187 = _tmp39 * _tmp60;
  const Scalar _tmp188 = _tmp40 * _tmp60;
  const Scalar _tmp189 = _tmp188 * _tmp86;
  const Scalar _tmp190 =
      -_tmp100 + _tmp102 * _tmp175 - _tmp103 - _tmp174 * _tmp29 + _tmp175 * _tmp34;
  const Scalar _tmp191 = _tmp104 * _tmp178 - _tmp190 * _tmp60;
  const Scalar _tmp192 = std::pow(_tmp105, Scalar(-2));
  const Scalar _tmp193 = _tmp191 * _tmp192;
  const Scalar _tmp194 = Scalar(1.0) * _tmp188;
  const Scalar _tmp195 = _tmp137 * _tmp40;
  const Scalar _tmp196 = Scalar(1.0) * _tmp55;
  const Scalar _tmp197 = _tmp132 * _tmp40;
  const Scalar _tmp198 = _tmp175 * _tmp59;
  const Scalar _tmp199 = _tmp131 * _tmp56;
  const Scalar _tmp200 = _tmp176 * _tmp54;
  const Scalar _tmp201 = _tmp124 * _tmp180 - _tmp200 * _tmp74;
  const Scalar _tmp202 = _tmp124 * _tmp182 - _tmp200 * _tmp73 - _tmp201 * _tmp82;
  const Scalar _tmp203 = _tmp105 * _tmp126;
  const Scalar _tmp204 = _tmp109 * (-_tmp104 * _tmp200 - _tmp106 * _tmp202 + _tmp124 * _tmp190 -
                                    _tmp127 * _tmp191 + _tmp185 * _tmp203);
  const Scalar _tmp205 = _tmp128 * _tmp83;
  const Scalar _tmp206 = _tmp193 * _tmp205;
  const Scalar _tmp207 = _tmp129 * _tmp183;
  const Scalar _tmp208 = _tmp202 + _tmp204 - _tmp206 + _tmp207;
  const Scalar _tmp209 = -_tmp112 * _tmp208 + _tmp185 * _tmp199;
  const Scalar _tmp210 = _tmp114 * _tmp31;
  const Scalar _tmp211 = _tmp113 * _tmp176;
  const Scalar _tmp212 = _tmp111 * _tmp56;
  const Scalar _tmp213 = _tmp107 * _tmp83;
  const Scalar _tmp214 = _tmp193 * _tmp213;
  const Scalar _tmp215 = _tmp84 * _tmp98;
  const Scalar _tmp216 = Scalar(1.0) * _tmp104;
  const Scalar _tmp217 = _tmp74 * _tmp96;
  const Scalar _tmp218 = Scalar(1.0) * _tmp73;
  const Scalar _tmp219 =
      -_tmp176 * _tmp217 + _tmp176 * _tmp218 + _tmp180 * _tmp97 - _tmp182 * _tmp94;
  const Scalar _tmp220 = _tmp105 * _tmp98;
  const Scalar _tmp221 = _tmp109 * (-_tmp106 * _tmp219 + _tmp176 * _tmp216 + _tmp185 * _tmp220 -
                                    _tmp190 * _tmp94 - _tmp191 * _tmp215);
  const Scalar _tmp222 = _tmp107 * _tmp108;
  const Scalar _tmp223 = _tmp183 * _tmp222;
  const Scalar _tmp224 = -_tmp214 + _tmp219 + _tmp221 + _tmp223;
  const Scalar _tmp225 = -_tmp112 * _tmp224 + _tmp185 * _tmp212;
  const Scalar _tmp226 = -_tmp123 * (-_tmp111 * _tmp186 + _tmp113 * _tmp198 - _tmp115 +
                                     _tmp210 * _tmp225 - _tmp211 * _tmp40 + _tmp224 * _tmp89) -
                         _tmp136 * (-_tmp131 * _tmp186 + _tmp132 * _tmp198 - _tmp133 -
                                    _tmp176 * _tmp197 + _tmp208 * _tmp89 + _tmp209 * _tmp210) -
                         _tmp140 * (_tmp137 * _tmp187 - _tmp138 * _tmp175 + _tmp178 * _tmp195 +
                                    _tmp193 * _tmp194 - _tmp193 * _tmp196) -
                         _tmp93 * (-_tmp175 * _tmp88 + _tmp176 * _tmp177 + _tmp185 * _tmp189 -
                                   _tmp186 * _tmp86 + _tmp187 * _tmp87);
  const Scalar _tmp227 = Scalar(9.6622558468725703) * _tmp226;
  const Scalar _tmp228 = _tmp146 * _tmp83;
  const Scalar _tmp229 = _tmp193 * _tmp228;
  const Scalar _tmp230 = _tmp108 * _tmp146;
  const Scalar _tmp231 = _tmp183 * _tmp230;
  const Scalar _tmp232 = Scalar(1.0) * _tmp76;
  const Scalar _tmp233 =
      -_tmp137 * _tmp181 - _tmp143 * _tmp229 + _tmp143 * _tmp231 + _tmp193 * _tmp232;
  const Scalar _tmp234 = _tmp185 * _tmp76;
  const Scalar _tmp235 = _tmp181 * _tmp84;
  const Scalar _tmp236 = _tmp131 * _tmp234 - _tmp131 * _tmp235 - _tmp150 * _tmp208 +
                         _tmp151 * _tmp204 - _tmp151 * _tmp206 + _tmp151 * _tmp207 + _tmp201;
  const Scalar _tmp237 = Scalar(1.0) * _tmp74;
  const Scalar _tmp238 = _tmp111 * _tmp234 - _tmp111 * _tmp235 - _tmp150 * _tmp224 -
                         _tmp151 * _tmp214 + _tmp151 * _tmp221 + _tmp151 * _tmp223 +
                         _tmp176 * _tmp237 - _tmp180 * _tmp94;
  const Scalar _tmp239 = _tmp234 * _tmp86 - _tmp235 * _tmp86;
  const Scalar _tmp240 = _tmp90 * _tmp95;
  const Scalar _tmp241 = std::pow(_tmp141, Scalar(-2));
  const Scalar _tmp242 = _tmp164 * _tmp241;
  const Scalar _tmp243 =
      std::pow(Scalar(std::pow(_tmp164, Scalar(2)) * _tmp241 + 1), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp244 = _tmp243 * (_tmp142 * (_tmp149 * (-_tmp229 + _tmp231 - _tmp233 * _tmp95) +
                                               _tmp153 * (_tmp146 * _tmp204 - _tmp146 * _tmp206 +
                                                          _tmp146 * _tmp207 - _tmp236 * _tmp95) +
                                               _tmp159 * (-_tmp146 * _tmp214 + _tmp146 * _tmp221 +
                                                          _tmp146 * _tmp223 - _tmp238 * _tmp95) -
                                               _tmp239 * _tmp240) -
                                    _tmp226 * _tmp242);
  const Scalar _tmp245 = Scalar(1.0) * std::cosh(_tmp170);
  const Scalar _tmp246 = Scalar(0.1034955) * _tmp167 * _tmp241;
  const Scalar _tmp247 = std::cosh(_tmp169);
  const Scalar _tmp248 = _tmp86 * _tmp90;
  const Scalar _tmp249 = _tmp248 * _tmp84;
  const Scalar _tmp250 = _tmp137 * _tmp139;
  const Scalar _tmp251 = _tmp132 * _tmp59;
  const Scalar _tmp252 = _tmp113 * _tmp59;
  const Scalar _tmp253 =
      _tmp122 * _tmp252 + _tmp135 * _tmp251 - _tmp249 * _tmp60 - _tmp250 * _tmp60;
  const Scalar _tmp254 = Scalar(1.0) / (_tmp253);
  const Scalar _tmp255 = _tmp139 * _tmp79;
  const Scalar _tmp256 = _tmp146 * _tmp163;
  const Scalar _tmp257 = _tmp79 * _tmp90;
  const Scalar _tmp258 = _tmp135 * _tmp79;
  const Scalar _tmp259 = _tmp122 * _tmp79;
  const Scalar _tmp260 = _tmp148 * _tmp255 + _tmp152 * _tmp258 + _tmp155 * _tmp157 * _tmp79 +
                         _tmp158 * _tmp259 + _tmp161 * _tmp257 - _tmp162 * _tmp256;
  const Scalar _tmp261 = std::asinh(_tmp254 * _tmp260);
  const Scalar _tmp262 = Scalar(1.0) * _tmp261;
  const Scalar _tmp263 = Scalar(9.6622558468725703) * _tmp253;
  const Scalar _tmp264 = -_tmp35 + p_c(1, 0);
  const Scalar _tmp265 = -_tmp30 + p_c(0, 0);
  const Scalar _tmp266 =
      std::sqrt(Scalar(std::pow(_tmp264, Scalar(2)) + std::pow(_tmp265, Scalar(2))));
  const Scalar _tmp267 = -_tmp261 * _tmp263 - _tmp266;
  const Scalar _tmp268 = Scalar(0.1034955) * _tmp254;
  const Scalar _tmp269 = _tmp267 * _tmp268;
  const Scalar _tmp270 = -std::sinh(_tmp262) - std::sinh(_tmp269);
  const Scalar _tmp271 = _tmp149 * _tmp193;
  const Scalar _tmp272 = _tmp112 * _tmp248;
  const Scalar _tmp273 = _tmp122 * _tmp59;
  const Scalar _tmp274 = _tmp132 * _tmp135;
  const Scalar _tmp275 = _tmp135 * _tmp59;
  const Scalar _tmp276 = _tmp185 * _tmp248;
  const Scalar _tmp277 = -_tmp122 * _tmp211 + _tmp176 * _tmp272 - _tmp176 * _tmp274 +
                         _tmp178 * _tmp250 + _tmp209 * _tmp275 + _tmp225 * _tmp273 +
                         _tmp271 * _tmp60 + _tmp276 * _tmp60;
  const Scalar _tmp278 = Scalar(9.6622558468725703) * _tmp277;
  const Scalar _tmp279 = std::cosh(_tmp269);
  const Scalar _tmp280 = Scalar(1.0) / (_tmp266);
  const Scalar _tmp281 = std::pow(_tmp253, Scalar(-2));
  const Scalar _tmp282 =
      std::pow(Scalar(std::pow(_tmp260, Scalar(2)) * _tmp281 + 1), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp283 = _tmp260 * _tmp281;
  const Scalar _tmp284 =
      _tmp282 *
      (_tmp254 * (_tmp233 * _tmp255 + _tmp236 * _tmp258 + _tmp238 * _tmp259 + _tmp239 * _tmp257) -
       _tmp277 * _tmp283);
  const Scalar _tmp285 = Scalar(0.1034955) * _tmp267 * _tmp281;
  const Scalar _tmp286 = Scalar(1.0) * std::cosh(_tmp262);
  const Scalar _tmp287 = _tmp122 * _tmp84;
  const Scalar _tmp288 = _tmp135 * _tmp84;
  const Scalar _tmp289 = _tmp111 * _tmp287 + _tmp131 * _tmp288 + _tmp249 + _tmp250;
  const Scalar _tmp290 = Scalar(1.0) / (_tmp289);
  const Scalar _tmp291 = _tmp122 * _tmp145;
  const Scalar _tmp292 = _tmp135 * _tmp145;
  const Scalar _tmp293 = -_tmp110 * _tmp291 - _tmp130 * _tmp292 - _tmp139 * _tmp147 -
                         _tmp155 * _tmp156 + _tmp160 * _tmp90 + _tmp256;
  const Scalar _tmp294 = std::asinh(_tmp290 * _tmp293);
  const Scalar _tmp295 = Scalar(9.6622558468725703) * _tmp289;
  const Scalar _tmp296 =
      -_tmp294 * _tmp295 - std::sqrt(Scalar(std::pow(Scalar(-_tmp43 + p_d(1, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp47 + p_d(0, 0)), Scalar(2))));
  const Scalar _tmp297 = Scalar(0.1034955) * _tmp290;
  const Scalar _tmp298 = _tmp296 * _tmp297;
  const Scalar _tmp299 = Scalar(1.0) * _tmp294;
  const Scalar _tmp300 = -std::sinh(_tmp298) - std::sinh(_tmp299);
  const Scalar _tmp301 = _tmp111 * _tmp122;
  const Scalar _tmp302 = _tmp131 * _tmp135;
  const Scalar _tmp303 = -_tmp185 * _tmp301 - _tmp185 * _tmp302 + _tmp208 * _tmp288 +
                         _tmp224 * _tmp287 - _tmp271 - _tmp276;
  const Scalar _tmp304 = Scalar(9.6622558468725703) * _tmp303;
  const Scalar _tmp305 = Scalar(1.0) * std::cosh(_tmp299);
  const Scalar _tmp306 = std::pow(_tmp289, Scalar(-2));
  const Scalar _tmp307 =
      std::pow(Scalar(std::pow(_tmp293, Scalar(2)) * _tmp306 + 1), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp308 = _tmp293 * _tmp306;
  const Scalar _tmp309 =
      _tmp307 *
      (_tmp290 * (_tmp139 * _tmp229 - _tmp139 * _tmp231 - _tmp204 * _tmp292 + _tmp206 * _tmp292 -
                  _tmp207 * _tmp292 + _tmp214 * _tmp291 - _tmp221 * _tmp291 - _tmp223 * _tmp291) -
       _tmp303 * _tmp308);
  const Scalar _tmp310 = Scalar(0.1034955) * _tmp296 * _tmp306;
  const Scalar _tmp311 = std::cosh(_tmp298);
  const Scalar _tmp312 = _tmp173 * _tmp37;
  const Scalar _tmp313 = -_tmp174 * _tmp75 + _tmp312 * _tmp69 - _tmp70;
  const Scalar _tmp314 = _tmp172 * (_tmp174 * _tmp54 - _tmp312 + _tmp39);
  const Scalar _tmp315 = _tmp314 * _tmp56;
  const Scalar _tmp316 = -_tmp313 * _tmp60 + _tmp315 * _tmp74;
  const Scalar _tmp317 = _tmp174 * _tmp68 - _tmp179;
  const Scalar _tmp318 = _tmp315 * _tmp73 - _tmp316 * _tmp82 - _tmp317 * _tmp60;
  const Scalar _tmp319 = _tmp184 * _tmp318;
  const Scalar _tmp320 = _tmp319 * _tmp55;
  const Scalar _tmp321 = _tmp102 * _tmp174 + _tmp174 * _tmp34 - _tmp29 * _tmp312 + _tmp99;
  const Scalar _tmp322 = _tmp104 * _tmp315 - _tmp321 * _tmp60;
  const Scalar _tmp323 = _tmp192 * _tmp322;
  const Scalar _tmp324 = _tmp213 * _tmp323;
  const Scalar _tmp325 =
      -_tmp217 * _tmp314 + _tmp218 * _tmp314 + _tmp313 * _tmp97 - _tmp317 * _tmp94;
  const Scalar _tmp326 = _tmp109 * (-_tmp106 * _tmp325 - _tmp215 * _tmp322 + _tmp216 * _tmp314 +
                                    _tmp220 * _tmp319 - _tmp321 * _tmp94);
  const Scalar _tmp327 = _tmp222 * _tmp318;
  const Scalar _tmp328 = -_tmp324 + _tmp325 + _tmp326 + _tmp327;
  const Scalar _tmp329 = -_tmp112 * _tmp328 + _tmp212 * _tmp319;
  const Scalar _tmp330 = _tmp113 * _tmp314;
  const Scalar _tmp331 = _tmp205 * _tmp323;
  const Scalar _tmp332 = _tmp314 * _tmp54;
  const Scalar _tmp333 = _tmp124 * _tmp313 - _tmp332 * _tmp74;
  const Scalar _tmp334 = _tmp124 * _tmp317 - _tmp332 * _tmp73 - _tmp333 * _tmp82;
  const Scalar _tmp335 = _tmp109 * (-_tmp104 * _tmp332 - _tmp106 * _tmp334 + _tmp124 * _tmp321 -
                                    _tmp127 * _tmp322 + _tmp203 * _tmp319);
  const Scalar _tmp336 = _tmp129 * _tmp318;
  const Scalar _tmp337 = -_tmp331 + _tmp334 + _tmp335 + _tmp336;
  const Scalar _tmp338 = -_tmp112 * _tmp337 + _tmp199 * _tmp319;
  const Scalar _tmp339 =
      -_tmp123 * (-_tmp111 * _tmp320 + _tmp174 * _tmp252 + _tmp210 * _tmp329 + _tmp328 * _tmp89 -
                  _tmp330 * _tmp40) -
      _tmp136 * (-_tmp131 * _tmp320 + _tmp174 * _tmp251 - _tmp197 * _tmp314 + _tmp210 * _tmp338 +
                 _tmp337 * _tmp89) -
      _tmp140 * (-_tmp138 * _tmp174 + _tmp194 * _tmp323 + _tmp195 * _tmp315 - _tmp196 * _tmp323) -
      _tmp93 * (-_tmp174 * _tmp88 + _tmp177 * _tmp314 + _tmp189 * _tmp319 - _tmp320 * _tmp86);
  const Scalar _tmp340 = _tmp228 * _tmp323;
  const Scalar _tmp341 = _tmp230 * _tmp318;
  const Scalar _tmp342 =
      -_tmp137 * _tmp316 - _tmp143 * _tmp340 + _tmp143 * _tmp341 + _tmp232 * _tmp323;
  const Scalar _tmp343 = _tmp316 * _tmp84;
  const Scalar _tmp344 = _tmp319 * _tmp76;
  const Scalar _tmp345 = -_tmp131 * _tmp343 + _tmp131 * _tmp344 - _tmp150 * _tmp337 -
                         _tmp151 * _tmp331 + _tmp151 * _tmp335 + _tmp151 * _tmp336 + _tmp333;
  const Scalar _tmp346 = -_tmp111 * _tmp343 + _tmp111 * _tmp344 - _tmp150 * _tmp328 -
                         _tmp151 * _tmp324 + _tmp151 * _tmp326 + _tmp151 * _tmp327 +
                         _tmp237 * _tmp314 - _tmp313 * _tmp94;
  const Scalar _tmp347 = -_tmp343 * _tmp86 + _tmp344 * _tmp86;
  const Scalar _tmp348 = _tmp243 * (_tmp142 * (_tmp149 * (-_tmp340 + _tmp341 - _tmp342 * _tmp95) +
                                               _tmp153 * (-_tmp146 * _tmp331 + _tmp146 * _tmp335 +
                                                          _tmp146 * _tmp336 - _tmp345 * _tmp95) +
                                               _tmp159 * (-_tmp146 * _tmp324 + _tmp146 * _tmp326 +
                                                          _tmp146 * _tmp327 - _tmp346 * _tmp95) -
                                               _tmp240 * _tmp347) -
                                    _tmp242 * _tmp339);
  const Scalar _tmp349 = Scalar(9.6622558468725703) * _tmp339;
  const Scalar _tmp350 = _tmp149 * _tmp323;
  const Scalar _tmp351 = _tmp248 * _tmp319;
  const Scalar _tmp352 = -_tmp122 * _tmp330 + _tmp250 * _tmp315 + _tmp272 * _tmp314 +
                         _tmp273 * _tmp329 - _tmp274 * _tmp314 + _tmp275 * _tmp338 +
                         _tmp350 * _tmp60 + _tmp351 * _tmp60;
  const Scalar _tmp353 = Scalar(9.6622558468725703) * _tmp352;
  const Scalar _tmp354 =
      _tmp282 *
      (_tmp254 * (_tmp255 * _tmp342 + _tmp257 * _tmp347 + _tmp258 * _tmp345 + _tmp259 * _tmp346) -
       _tmp283 * _tmp352);
  const Scalar _tmp355 = _tmp287 * _tmp328 + _tmp288 * _tmp337 - _tmp301 * _tmp319 -
                         _tmp302 * _tmp319 - _tmp350 - _tmp351;
  const Scalar _tmp356 = Scalar(9.6622558468725703) * _tmp355;
  const Scalar _tmp357 =
      _tmp307 *
      (_tmp290 * (_tmp139 * _tmp340 - _tmp139 * _tmp341 + _tmp291 * _tmp324 - _tmp291 * _tmp326 -
                  _tmp291 * _tmp327 + _tmp292 * _tmp331 - _tmp292 * _tmp335 - _tmp292 * _tmp336) -
       _tmp308 * _tmp355);

  // Output terms (1)
  Eigen::Matrix<Scalar, 4, 3> _res;

  _res(0, 0) = 0;
  _res(1, 0) =
      _tmp166 *
          (-_tmp244 * _tmp245 -
           _tmp247 * (_tmp168 * (-_tmp165 * _tmp227 - _tmp166 * _tmp244) - _tmp226 * _tmp246)) +
      _tmp171 * _tmp227;
  _res(2, 0) =
      _tmp263 *
          (-_tmp279 * (_tmp268 * (-_tmp261 * _tmp278 - _tmp263 * _tmp284 - _tmp265 * _tmp280) -
                       _tmp277 * _tmp285) -
           _tmp284 * _tmp286) +
      _tmp270 * _tmp278;
  _res(3, 0) =
      _tmp295 *
          (-_tmp305 * _tmp309 -
           _tmp311 * (_tmp297 * (-_tmp294 * _tmp304 - _tmp295 * _tmp309) - _tmp303 * _tmp310)) +
      _tmp300 * _tmp304;
  _res(0, 1) = 0;
  _res(1, 1) =
      _tmp166 *
          (-_tmp245 * _tmp348 -
           _tmp247 * (_tmp168 * (-_tmp165 * _tmp349 - _tmp166 * _tmp348) - _tmp246 * _tmp339)) +
      _tmp171 * _tmp349;
  _res(2, 1) =
      _tmp263 *
          (-_tmp279 * (_tmp268 * (-_tmp261 * _tmp353 - _tmp263 * _tmp354 - _tmp264 * _tmp280) -
                       _tmp285 * _tmp352) -
           _tmp286 * _tmp354) +
      _tmp270 * _tmp353;
  _res(3, 1) =
      _tmp295 *
          (-_tmp305 * _tmp357 -
           _tmp311 * (_tmp297 * (-_tmp294 * _tmp356 - _tmp295 * _tmp357) - _tmp310 * _tmp355)) +
      _tmp300 * _tmp356;
  _res(0, 2) = 0;
  _res(1, 2) = 0;
  _res(2, 2) = 0;
  _res(3, 2) = 0;

  return _res;
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym

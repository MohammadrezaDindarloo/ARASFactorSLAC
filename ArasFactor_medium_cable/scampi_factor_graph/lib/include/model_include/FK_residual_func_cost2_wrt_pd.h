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
 * Symbolic function: FK_residual_func_cost2_wrt_pd
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
Eigen::Matrix<Scalar, 4, 3> FkResidualFuncCost2WrtPd(
    const Scalar fh1, const Scalar fv1, const sym::Rot3<Scalar>& DeltaRot,
    const sym::Pose3<Scalar>& TransformationMatrix, const Eigen::Matrix<Scalar, 4, 1>& encoder,
    const Eigen::Matrix<Scalar, 3, 1>& p_a, const Eigen::Matrix<Scalar, 3, 1>& p_b,
    const Eigen::Matrix<Scalar, 3, 1>& p_c, const Eigen::Matrix<Scalar, 3, 1>& p_d,
    const Scalar epsilon) {
  // Total ops: 1045

  // Unused inputs
  (void)encoder;
  (void)epsilon;

  // Input arrays
  const Eigen::Matrix<Scalar, 4, 1>& _DeltaRot = DeltaRot.Data();
  const Eigen::Matrix<Scalar, 7, 1>& _TransformationMatrix = TransformationMatrix.Data();

  // Intermediate terms (321)
  const Scalar _tmp0 =
      _DeltaRot[0] * _TransformationMatrix[2] + _DeltaRot[1] * _TransformationMatrix[3] -
      _DeltaRot[2] * _TransformationMatrix[0] + _DeltaRot[3] * _TransformationMatrix[1];
  const Scalar _tmp1 =
      _DeltaRot[0] * _TransformationMatrix[3] - _DeltaRot[1] * _TransformationMatrix[2] +
      _DeltaRot[2] * _TransformationMatrix[1] + _DeltaRot[3] * _TransformationMatrix[0];
  const Scalar _tmp2 = 2 * _tmp1;
  const Scalar _tmp3 = _tmp0 * _tmp2;
  const Scalar _tmp4 =
      -_DeltaRot[0] * _TransformationMatrix[1] + _DeltaRot[1] * _TransformationMatrix[0] +
      _DeltaRot[2] * _TransformationMatrix[3] + _DeltaRot[3] * _TransformationMatrix[2];
  const Scalar _tmp5 =
      -2 * _DeltaRot[0] * _TransformationMatrix[0] - 2 * _DeltaRot[1] * _TransformationMatrix[1] -
      2 * _DeltaRot[2] * _TransformationMatrix[2] + 2 * _DeltaRot[3] * _TransformationMatrix[3];
  const Scalar _tmp6 = _tmp4 * _tmp5;
  const Scalar _tmp7 = Scalar(0.20999999999999999) * _tmp3 + Scalar(0.20999999999999999) * _tmp6;
  const Scalar _tmp8 = 2 * _tmp0 * _tmp4;
  const Scalar _tmp9 = _tmp1 * _tmp5;
  const Scalar _tmp10 = _tmp8 - _tmp9;
  const Scalar _tmp11 = -Scalar(0.010999999999999999) * _tmp10;
  const Scalar _tmp12 = -2 * std::pow(_tmp1, Scalar(2));
  const Scalar _tmp13 = -2 * std::pow(_tmp4, Scalar(2));
  const Scalar _tmp14 = Scalar(0.20999999999999999) * _tmp12 +
                        Scalar(0.20999999999999999) * _tmp13 + Scalar(0.20999999999999999);
  const Scalar _tmp15 = _tmp11 - _tmp14;
  const Scalar _tmp16 = _tmp15 + _tmp7;
  const Scalar _tmp17 = _TransformationMatrix[5] + _tmp16;
  const Scalar _tmp18 = _tmp17 - p_b(1, 0);
  const Scalar _tmp19 = 1 - 2 * std::pow(_tmp0, Scalar(2));
  const Scalar _tmp20 = Scalar(0.20999999999999999) * _tmp13 + Scalar(0.20999999999999999) * _tmp19;
  const Scalar _tmp21 = _tmp2 * _tmp4;
  const Scalar _tmp22 = _tmp0 * _tmp5;
  const Scalar _tmp23 = _tmp21 + _tmp22;
  const Scalar _tmp24 = -Scalar(0.010999999999999999) * _tmp23;
  const Scalar _tmp25 = Scalar(0.20999999999999999) * _tmp3 - Scalar(0.20999999999999999) * _tmp6;
  const Scalar _tmp26 = _tmp24 - _tmp25;
  const Scalar _tmp27 = _tmp20 + _tmp26;
  const Scalar _tmp28 = _TransformationMatrix[4] + _tmp27;
  const Scalar _tmp29 = _tmp28 - p_b(0, 0);
  const Scalar _tmp30 = Scalar(1.0) / (_tmp29);
  const Scalar _tmp31 = _tmp18 * _tmp30;
  const Scalar _tmp32 = _tmp24 + _tmp25;
  const Scalar _tmp33 = _tmp20 + _tmp32;
  const Scalar _tmp34 = _TransformationMatrix[4] + _tmp33;
  const Scalar _tmp35 = _tmp34 - p_c(0, 0);
  const Scalar _tmp36 = _tmp11 + _tmp14;
  const Scalar _tmp37 = _tmp36 + _tmp7;
  const Scalar _tmp38 = _TransformationMatrix[5] + _tmp37;
  const Scalar _tmp39 = _tmp38 - p_c(1, 0);
  const Scalar _tmp40 = std::pow(Scalar(std::pow(_tmp35, Scalar(2)) + std::pow(_tmp39, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp41 = _tmp35 * _tmp40;
  const Scalar _tmp42 = _tmp39 * _tmp40;
  const Scalar _tmp43 = Scalar(1.0) / (_tmp31 * _tmp41 - _tmp42);
  const Scalar _tmp44 = _tmp41 * _tmp43;
  const Scalar _tmp45 = Scalar(1.0) * _tmp27;
  const Scalar _tmp46 = Scalar(1.0) * _tmp16;
  const Scalar _tmp47 = -_tmp46;
  const Scalar _tmp48 = Scalar(1.0) / (_tmp37 + _tmp47);
  const Scalar _tmp49 = -_tmp33 + _tmp45;
  const Scalar _tmp50 = _tmp48 * _tmp49;
  const Scalar _tmp51 = _tmp45 + _tmp46 * _tmp50;
  const Scalar _tmp52 = 0;
  const Scalar _tmp53 = -_tmp20;
  const Scalar _tmp54 = _tmp32 + _tmp53;
  const Scalar _tmp55 = _TransformationMatrix[4] + _tmp54;
  const Scalar _tmp56 = _tmp55 - p_d(0, 0);
  const Scalar _tmp57 = Scalar(0.20999999999999999) * _tmp21 - Scalar(0.20999999999999999) * _tmp22;
  const Scalar _tmp58 = -_tmp57;
  const Scalar _tmp59 =
      -Scalar(0.010999999999999999) * _tmp12 - Scalar(0.010999999999999999) * _tmp19;
  const Scalar _tmp60 = Scalar(0.20999999999999999) * _tmp8 + Scalar(0.20999999999999999) * _tmp9;
  const Scalar _tmp61 = _tmp59 + _tmp60;
  const Scalar _tmp62 = _tmp58 + _tmp61;
  const Scalar _tmp63 = -_tmp7;
  const Scalar _tmp64 = _tmp36 + _tmp63;
  const Scalar _tmp65 = _TransformationMatrix[5] + _tmp64;
  const Scalar _tmp66 = _tmp65 - p_d(1, 0);
  const Scalar _tmp67 = std::pow(_tmp66, Scalar(2));
  const Scalar _tmp68 = std::pow(_tmp56, Scalar(2));
  const Scalar _tmp69 = _tmp67 + _tmp68;
  const Scalar _tmp70 = std::pow(_tmp69, Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp71 = _tmp62 * _tmp70;
  const Scalar _tmp72 = _tmp59 - _tmp60;
  const Scalar _tmp73 = _tmp57 + _tmp72;
  const Scalar _tmp74 = _tmp70 * _tmp73;
  const Scalar _tmp75 = _tmp56 * _tmp74;
  const Scalar _tmp76 = _tmp66 * _tmp70;
  const Scalar _tmp77 = _tmp31 * _tmp70;
  const Scalar _tmp78 = _tmp56 * _tmp77 - _tmp76;
  const Scalar _tmp79 = _tmp57 + _tmp61;
  const Scalar _tmp80 = _tmp43 * (_tmp41 * _tmp73 - _tmp41 * _tmp79);
  const Scalar _tmp81 = _tmp31 * _tmp73;
  const Scalar _tmp82 = _tmp43 * (-_tmp41 * _tmp81 + _tmp42 * _tmp79);
  const Scalar _tmp83 = -_tmp31 * _tmp75 + _tmp62 * _tmp76 - _tmp78 * _tmp82;
  const Scalar _tmp84 = -_tmp50 * _tmp83 - _tmp56 * _tmp71 + _tmp75 - _tmp78 * _tmp80;
  const Scalar _tmp85 = Scalar(1.0) / (_tmp84);
  const Scalar _tmp86 = _tmp78 * _tmp85;
  const Scalar _tmp87 = _tmp52 * _tmp86;
  const Scalar _tmp88 = _tmp52 * _tmp85;
  const Scalar _tmp89 = _tmp70 * _tmp88;
  const Scalar _tmp90 = Scalar(43.164000000000001) - fv1;
  const Scalar _tmp91 =
      std::sqrt(Scalar(std::pow(_tmp18, Scalar(2)) + std::pow(_tmp29, Scalar(2))));
  const Scalar _tmp92 = _tmp30 * _tmp91;
  const Scalar _tmp93 = _tmp90 * _tmp92;
  const Scalar _tmp94 = Scalar(1.0) * _tmp48;
  const Scalar _tmp95 = _tmp49 * _tmp82 * _tmp94 - Scalar(1.0) * _tmp80;
  const Scalar _tmp96 = _tmp64 * _tmp70;
  const Scalar _tmp97 = Scalar(1.0) / (_tmp91);
  const Scalar _tmp98 = _tmp92 * (-_tmp16 * _tmp29 * _tmp97 + _tmp18 * _tmp27 * _tmp97);
  const Scalar _tmp99 = _tmp43 * (-_tmp33 * _tmp42 + _tmp37 * _tmp41 + _tmp41 * _tmp98);
  const Scalar _tmp100 = _tmp56 * _tmp70;
  const Scalar _tmp101 = _tmp100 * _tmp98 - _tmp54 * _tmp76 + _tmp56 * _tmp96 - _tmp78 * _tmp99;
  const Scalar _tmp102 = _tmp101 * _tmp85;
  const Scalar _tmp103 = -_tmp102 * _tmp95 - Scalar(1.0) * _tmp99;
  const Scalar _tmp104 = Scalar(1.0) / (_tmp101);
  const Scalar _tmp105 = _tmp103 * _tmp104;
  const Scalar _tmp106 = _tmp105 * _tmp84;
  const Scalar _tmp107 = _tmp106 + _tmp95;
  const Scalar _tmp108 = -_tmp107 * _tmp86 + Scalar(1.0);
  const Scalar _tmp109 = _tmp107 * _tmp85;
  const Scalar _tmp110 = _tmp109 * _tmp70;
  const Scalar _tmp111 = _tmp26 + _tmp53;
  const Scalar _tmp112 = _TransformationMatrix[4] + _tmp111 - p_a(0, 0);
  const Scalar _tmp113 = _tmp15 + _tmp63;
  const Scalar _tmp114 = _TransformationMatrix[5] + _tmp113 - p_a(1, 0);
  const Scalar _tmp115 =
      std::pow(Scalar(std::pow(_tmp112, Scalar(2)) + std::pow(_tmp114, Scalar(2))),
               Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp116 = _tmp114 * _tmp115;
  const Scalar _tmp117 = _tmp116 * fh1;
  const Scalar _tmp118 = _tmp117 * _tmp92;
  const Scalar _tmp119 = _tmp31 * _tmp82 + _tmp81;
  const Scalar _tmp120 = -_tmp119 * _tmp50 + _tmp31 * _tmp80 - _tmp73;
  const Scalar _tmp121 = -_tmp102 * _tmp120 + _tmp31 * _tmp99 - _tmp98;
  const Scalar _tmp122 = _tmp104 * _tmp84;
  const Scalar _tmp123 = _tmp121 * _tmp122;
  const Scalar _tmp124 = _tmp120 + _tmp123;
  const Scalar _tmp125 = -_tmp124 * _tmp86 - _tmp31;
  const Scalar _tmp126 = _tmp124 * _tmp85;
  const Scalar _tmp127 = _tmp126 * _tmp70;
  const Scalar _tmp128 = _tmp112 * _tmp115;
  const Scalar _tmp129 = _tmp128 * fh1;
  const Scalar _tmp130 = _tmp129 * _tmp92;
  const Scalar _tmp131 = Scalar(1.0) * _tmp104;
  const Scalar _tmp132 = _tmp131 * _tmp70;
  const Scalar _tmp133 = _tmp131 * _tmp44;
  const Scalar _tmp134 = fh1 * (_tmp111 * _tmp116 - _tmp113 * _tmp128);
  const Scalar _tmp135 = _tmp134 * _tmp92;
  const Scalar _tmp136 = -_tmp118 * (_tmp108 * _tmp44 + _tmp110 * _tmp56) -
                         _tmp130 * (_tmp125 * _tmp44 + _tmp127 * _tmp56 + Scalar(1.0)) -
                         _tmp135 * (_tmp132 * _tmp56 - _tmp133 * _tmp78) -
                         _tmp93 * (-_tmp44 * _tmp87 + _tmp56 * _tmp89);
  const Scalar _tmp137 = Scalar(1.0) / (_tmp136);
  const Scalar _tmp138 = _tmp47 + _tmp64;
  const Scalar _tmp139 = _tmp138 * _tmp50;
  const Scalar _tmp140 = Scalar(1.0) / (-_tmp139 + _tmp45 - _tmp54);
  const Scalar _tmp141 = Scalar(1.0) * _tmp140;
  const Scalar _tmp142 = _tmp122 * _tmp141;
  const Scalar _tmp143 = -_tmp131 * _tmp83 + _tmp138 * _tmp142;
  const Scalar _tmp144 = Scalar(1.0) * _tmp134;
  const Scalar _tmp145 = _tmp83 * _tmp85;
  const Scalar _tmp146 = _tmp138 * _tmp140;
  const Scalar _tmp147 = _tmp119 + _tmp123 * _tmp146 - _tmp124 * _tmp145;
  const Scalar _tmp148 = Scalar(1.0) * _tmp129;
  const Scalar _tmp149 = fh1 * (_tmp58 + _tmp72);
  const Scalar _tmp150 = -Scalar(5.1796800000000003) * _tmp10 - _tmp113 * fv1 - _tmp116 * _tmp149;
  const Scalar _tmp151 = _tmp141 * _tmp50;
  const Scalar _tmp152 = _tmp139 * _tmp141 + Scalar(1.0);
  const Scalar _tmp153 = _tmp106 * _tmp146 - _tmp107 * _tmp145 - Scalar(1.0) * _tmp82;
  const Scalar _tmp154 = Scalar(1.0) * _tmp117;
  const Scalar _tmp155 = _tmp140 * _tmp51;
  const Scalar _tmp156 = -_tmp138 * _tmp155 - _tmp145 * _tmp52 + _tmp47;
  const Scalar _tmp157 = _tmp138 * _tmp48;
  const Scalar _tmp158 = _tmp111 * fv1 + _tmp128 * _tmp149 + Scalar(5.1796800000000003) * _tmp23;
  const Scalar _tmp159 =
      _tmp144 * (_tmp142 - _tmp143 * _tmp94) + _tmp148 * (_tmp123 * _tmp141 - _tmp147 * _tmp94) +
      Scalar(1.0) * _tmp150 * (_tmp151 - _tmp152 * _tmp94) +
      _tmp154 * (_tmp106 * _tmp141 - _tmp153 * _tmp94) +
      Scalar(1.0) * _tmp158 * (_tmp141 * _tmp157 - _tmp141) +
      Scalar(1.0) * _tmp90 * (-_tmp141 * _tmp51 - _tmp156 * _tmp94 + Scalar(1.0));
  const Scalar _tmp160 = std::asinh(_tmp137 * _tmp159);
  const Scalar _tmp161 = Scalar(1.0) * _tmp160;
  const Scalar _tmp162 = Scalar(1.0) * std::cosh(_tmp161);
  const Scalar _tmp163 = std::pow(_tmp136, Scalar(-2));
  const Scalar _tmp164 =
      std::pow(Scalar(std::pow(_tmp159, Scalar(2)) * _tmp163 + 1), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp165 = std::pow(_tmp69, Scalar(Scalar(-3) / Scalar(2)));
  const Scalar _tmp166 = _tmp165 * _tmp56 * _tmp66;
  const Scalar _tmp167 = _tmp165 * _tmp68;
  const Scalar _tmp168 = -_tmp166 + _tmp167 * _tmp31 - _tmp77;
  const Scalar _tmp169 = -_tmp166 * _tmp54 + _tmp167 * _tmp64 + _tmp167 * _tmp98 -
                         _tmp168 * _tmp99 - _tmp70 * _tmp98 - _tmp96;
  const Scalar _tmp170 = _tmp169 * _tmp85;
  const Scalar _tmp171 = _tmp166 * _tmp62;
  const Scalar _tmp172 = -_tmp167 * _tmp81 - _tmp168 * _tmp82 + _tmp171 + _tmp31 * _tmp74;
  const Scalar _tmp173 =
      -_tmp167 * _tmp62 + _tmp167 * _tmp73 - _tmp168 * _tmp80 - _tmp172 * _tmp50 + _tmp71 - _tmp74;
  const Scalar _tmp174 = std::pow(_tmp84, Scalar(-2));
  const Scalar _tmp175 = _tmp173 * _tmp174;
  const Scalar _tmp176 = _tmp101 * _tmp175;
  const Scalar _tmp177 = _tmp122 * (-_tmp170 * _tmp95 + _tmp176 * _tmp95);
  const Scalar _tmp178 = std::pow(_tmp101, Scalar(-2));
  const Scalar _tmp179 = _tmp169 * _tmp178;
  const Scalar _tmp180 = _tmp179 * _tmp84;
  const Scalar _tmp181 = _tmp103 * _tmp180;
  const Scalar _tmp182 = _tmp105 * _tmp173;
  const Scalar _tmp183 = _tmp177 - _tmp181 + _tmp182;
  const Scalar _tmp184 = _tmp183 * _tmp85;
  const Scalar _tmp185 = _tmp107 * _tmp175;
  const Scalar _tmp186 = -_tmp109 * _tmp168 - _tmp183 * _tmp86 + _tmp185 * _tmp78;
  const Scalar _tmp187 = _tmp44 * _tmp78;
  const Scalar _tmp188 = Scalar(1.0) * _tmp187;
  const Scalar _tmp189 = Scalar(1.0) * _tmp100;
  const Scalar _tmp190 = _tmp168 * _tmp44;
  const Scalar _tmp191 = _tmp104 * _tmp121;
  const Scalar _tmp192 = _tmp173 * _tmp191;
  const Scalar _tmp193 = _tmp122 * (-_tmp120 * _tmp170 + _tmp120 * _tmp176);
  const Scalar _tmp194 = _tmp121 * _tmp180;
  const Scalar _tmp195 = _tmp192 + _tmp193 - _tmp194;
  const Scalar _tmp196 = _tmp100 * _tmp85;
  const Scalar _tmp197 = _tmp124 * _tmp175;
  const Scalar _tmp198 = -_tmp126 * _tmp168 - _tmp195 * _tmp86 + _tmp197 * _tmp78;
  const Scalar _tmp199 = _tmp175 * _tmp52;
  const Scalar _tmp200 =
      -_tmp118 *
          (_tmp100 * _tmp184 - _tmp100 * _tmp185 + _tmp109 * _tmp167 - _tmp110 + _tmp186 * _tmp44) -
      _tmp130 * (-_tmp100 * _tmp197 + _tmp126 * _tmp167 - _tmp127 + _tmp195 * _tmp196 +
                 _tmp198 * _tmp44) -
      _tmp135 * (_tmp131 * _tmp167 - _tmp131 * _tmp190 - _tmp132 + _tmp179 * _tmp188 -
                 _tmp179 * _tmp189) -
      _tmp93 *
          (-_tmp100 * _tmp199 + _tmp167 * _tmp88 + _tmp187 * _tmp199 - _tmp190 * _tmp88 - _tmp89);
  const Scalar _tmp201 = _tmp159 * _tmp163;
  const Scalar _tmp202 = -_tmp172 * _tmp88 + _tmp199 * _tmp83;
  const Scalar _tmp203 = _tmp90 * _tmp94;
  const Scalar _tmp204 = -_tmp109 * _tmp172 - _tmp145 * _tmp183 + _tmp146 * _tmp177 -
                         _tmp146 * _tmp181 + _tmp146 * _tmp182 + _tmp185 * _tmp83;
  const Scalar _tmp205 = -_tmp126 * _tmp172 - _tmp145 * _tmp195 + _tmp146 * _tmp192 +
                         _tmp146 * _tmp193 - _tmp146 * _tmp194 + _tmp197 * _tmp83;
  const Scalar _tmp206 = _tmp141 * _tmp180;
  const Scalar _tmp207 = Scalar(1.0) * _tmp83;
  const Scalar _tmp208 = _tmp104 * _tmp141;
  const Scalar _tmp209 = _tmp173 * _tmp208;
  const Scalar _tmp210 =
      -_tmp131 * _tmp172 - _tmp138 * _tmp206 + _tmp138 * _tmp209 + _tmp179 * _tmp207;
  const Scalar _tmp211 = _tmp164 * (_tmp137 * (_tmp144 * (-_tmp206 + _tmp209 - _tmp210 * _tmp94) +
                                               _tmp148 * (_tmp141 * _tmp192 + _tmp141 * _tmp193 -
                                                          _tmp141 * _tmp194 - _tmp205 * _tmp94) +
                                               _tmp154 * (_tmp141 * _tmp177 - _tmp141 * _tmp181 +
                                                          _tmp141 * _tmp182 - _tmp204 * _tmp94) -
                                               _tmp202 * _tmp203) -
                                    _tmp200 * _tmp201);
  const Scalar _tmp212 = Scalar(9.6622558468725703) * _tmp160;
  const Scalar _tmp213 = Scalar(9.6622558468725703) * _tmp136;
  const Scalar _tmp214 = Scalar(0.1034955) * _tmp137;
  const Scalar _tmp215 =
      -_tmp160 * _tmp213 - std::sqrt(Scalar(std::pow(Scalar(-_tmp17 + p_b(1, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp28 + p_b(0, 0)), Scalar(2))));
  const Scalar _tmp216 = Scalar(0.1034955) * _tmp163 * _tmp215;
  const Scalar _tmp217 = _tmp214 * _tmp215;
  const Scalar _tmp218 = std::cosh(_tmp217);
  const Scalar _tmp219 = -Scalar(9.6622558468725703) * std::sinh(_tmp161) -
                         Scalar(9.6622558468725703) * std::sinh(_tmp217);
  const Scalar _tmp220 = _tmp131 * _tmp134;
  const Scalar _tmp221 = _tmp220 * _tmp43;
  const Scalar _tmp222 = _tmp129 * _tmp43;
  const Scalar _tmp223 = _tmp117 * _tmp43;
  const Scalar _tmp224 =
      _tmp108 * _tmp223 + _tmp125 * _tmp222 - _tmp221 * _tmp78 - _tmp43 * _tmp87 * _tmp90;
  const Scalar _tmp225 = Scalar(1.0) / (_tmp224);
  const Scalar _tmp226 = _tmp134 * _tmp48;
  const Scalar _tmp227 = _tmp141 * _tmp158;
  const Scalar _tmp228 = _tmp48 * _tmp90;
  const Scalar _tmp229 = _tmp129 * _tmp48;
  const Scalar _tmp230 = _tmp117 * _tmp48;
  const Scalar _tmp231 = _tmp143 * _tmp226 + _tmp147 * _tmp229 + _tmp150 * _tmp152 * _tmp48 +
                         _tmp153 * _tmp230 + _tmp156 * _tmp228 - _tmp157 * _tmp227;
  const Scalar _tmp232 = std::asinh(_tmp225 * _tmp231);
  const Scalar _tmp233 = Scalar(1.0) * _tmp232;
  const Scalar _tmp234 = Scalar(9.6622558468725703) * _tmp224;
  const Scalar _tmp235 =
      -_tmp232 * _tmp234 - std::sqrt(Scalar(std::pow(Scalar(-_tmp34 + p_c(0, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp38 + p_c(1, 0)), Scalar(2))));
  const Scalar _tmp236 = Scalar(0.1034955) * _tmp225;
  const Scalar _tmp237 = _tmp235 * _tmp236;
  const Scalar _tmp238 = -std::sinh(_tmp233) - std::sinh(_tmp237);
  const Scalar _tmp239 = _tmp88 * _tmp90;
  const Scalar _tmp240 = _tmp239 * _tmp43;
  const Scalar _tmp241 = _tmp199 * _tmp90;
  const Scalar _tmp242 = _tmp43 * _tmp78;
  const Scalar _tmp243 = _tmp144 * _tmp179;
  const Scalar _tmp244 = -_tmp168 * _tmp221 - _tmp168 * _tmp240 + _tmp186 * _tmp223 +
                         _tmp198 * _tmp222 + _tmp241 * _tmp242 + _tmp242 * _tmp243;
  const Scalar _tmp245 = Scalar(9.6622558468725703) * _tmp244;
  const Scalar _tmp246 = std::pow(_tmp224, Scalar(-2));
  const Scalar _tmp247 =
      std::pow(Scalar(std::pow(_tmp231, Scalar(2)) * _tmp246 + 1), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp248 = _tmp231 * _tmp246;
  const Scalar _tmp249 =
      _tmp247 *
      (_tmp225 * (_tmp202 * _tmp228 + _tmp204 * _tmp230 + _tmp205 * _tmp229 + _tmp210 * _tmp226) -
       _tmp244 * _tmp248);
  const Scalar _tmp250 = Scalar(0.1034955) * _tmp235 * _tmp246;
  const Scalar _tmp251 = std::cosh(_tmp237);
  const Scalar _tmp252 = Scalar(1.0) * std::cosh(_tmp233);
  const Scalar _tmp253 = -_TransformationMatrix[6] - _tmp62 + p_d(2, 0);
  const Scalar _tmp254 = -_tmp55 + p_d(0, 0);
  const Scalar _tmp255 = -_tmp65 + p_d(1, 0);
  const Scalar _tmp256 = std::pow(_tmp254, Scalar(2)) + std::pow(_tmp255, Scalar(2));
  const Scalar _tmp257 =
      std::pow(Scalar(std::pow(_tmp253, Scalar(2)) + _tmp256), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp258 = _tmp109 * _tmp117 + _tmp126 * _tmp129 + _tmp220 + _tmp239;
  const Scalar _tmp259 = Scalar(1.0) / (_tmp258);
  const Scalar _tmp260 = _tmp117 * _tmp140;
  const Scalar _tmp261 = _tmp129 * _tmp140;
  const Scalar _tmp262 = -_tmp106 * _tmp260 - _tmp123 * _tmp261 - _tmp134 * _tmp142 -
                         _tmp150 * _tmp151 + _tmp155 * _tmp90 + _tmp227;
  const Scalar _tmp263 = std::asinh(_tmp259 * _tmp262);
  const Scalar _tmp264 = Scalar(1.0) * _tmp263;
  const Scalar _tmp265 = Scalar(1.0) * std::cosh(_tmp264);
  const Scalar _tmp266 = std::pow(_tmp258, Scalar(-2));
  const Scalar _tmp267 =
      std::pow(Scalar(std::pow(_tmp262, Scalar(2)) * _tmp266 + 1), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp268 = _tmp129 * _tmp85;
  const Scalar _tmp269 = _tmp117 * _tmp184 - _tmp117 * _tmp185 - _tmp129 * _tmp197 +
                         _tmp195 * _tmp268 - _tmp241 - _tmp243;
  const Scalar _tmp270 = _tmp262 * _tmp266;
  const Scalar _tmp271 =
      _tmp267 *
      (_tmp259 * (_tmp134 * _tmp206 - _tmp134 * _tmp209 - _tmp177 * _tmp260 + _tmp181 * _tmp260 -
                  _tmp182 * _tmp260 - _tmp192 * _tmp261 - _tmp193 * _tmp261 + _tmp194 * _tmp261) -
       _tmp269 * _tmp270);
  const Scalar _tmp272 = std::sqrt(_tmp256);
  const Scalar _tmp273 = Scalar(9.6622558468725703) * _tmp258;
  const Scalar _tmp274 = -_tmp263 * _tmp273 - _tmp272;
  const Scalar _tmp275 = Scalar(0.1034955) * _tmp259;
  const Scalar _tmp276 = _tmp274 * _tmp275;
  const Scalar _tmp277 = std::cosh(_tmp276);
  const Scalar _tmp278 = Scalar(0.1034955) * _tmp266 * _tmp274;
  const Scalar _tmp279 = Scalar(1.0) / (_tmp272);
  const Scalar _tmp280 = Scalar(9.6622558468725703) * _tmp269;
  const Scalar _tmp281 = -std::sinh(_tmp264) - std::sinh(_tmp276);
  const Scalar _tmp282 = _tmp165 * _tmp67;
  const Scalar _tmp283 = _tmp166 * _tmp31 - _tmp282 + _tmp70;
  const Scalar _tmp284 =
      _tmp166 * _tmp64 + _tmp166 * _tmp98 - _tmp282 * _tmp54 - _tmp283 * _tmp99 + _tmp54 * _tmp70;
  const Scalar _tmp285 = _tmp178 * _tmp284;
  const Scalar _tmp286 = _tmp284 * _tmp85;
  const Scalar _tmp287 = -_tmp166 * _tmp81 + _tmp282 * _tmp62 - _tmp283 * _tmp82 - _tmp71;
  const Scalar _tmp288 = _tmp166 * _tmp73 - _tmp171 - _tmp283 * _tmp80 - _tmp287 * _tmp50;
  const Scalar _tmp289 = _tmp174 * _tmp288;
  const Scalar _tmp290 = _tmp101 * _tmp289;
  const Scalar _tmp291 = _tmp122 * (-_tmp120 * _tmp286 + _tmp120 * _tmp290);
  const Scalar _tmp292 = _tmp191 * _tmp288;
  const Scalar _tmp293 = _tmp285 * _tmp84;
  const Scalar _tmp294 = _tmp121 * _tmp293;
  const Scalar _tmp295 = _tmp291 + _tmp292 - _tmp294;
  const Scalar _tmp296 = _tmp124 * _tmp289;
  const Scalar _tmp297 = -_tmp126 * _tmp283 - _tmp295 * _tmp86 + _tmp296 * _tmp78;
  const Scalar _tmp298 = _tmp289 * _tmp52;
  const Scalar _tmp299 = _tmp107 * _tmp289;
  const Scalar _tmp300 = _tmp103 * _tmp293;
  const Scalar _tmp301 = _tmp122 * (-_tmp286 * _tmp95 + _tmp290 * _tmp95);
  const Scalar _tmp302 = _tmp105 * _tmp288;
  const Scalar _tmp303 = -_tmp300 + _tmp301 + _tmp302;
  const Scalar _tmp304 = -_tmp109 * _tmp283 + _tmp299 * _tmp78 - _tmp303 * _tmp86;
  const Scalar _tmp305 =
      -_tmp118 * (-_tmp100 * _tmp299 + _tmp109 * _tmp166 + _tmp196 * _tmp303 + _tmp304 * _tmp44) -
      _tmp130 * (-_tmp100 * _tmp296 + _tmp126 * _tmp166 + _tmp196 * _tmp295 + _tmp297 * _tmp44) -
      _tmp135 * (_tmp131 * _tmp166 - _tmp133 * _tmp283 + _tmp188 * _tmp285 - _tmp189 * _tmp285) -
      _tmp93 *
          (-_tmp100 * _tmp298 + _tmp166 * _tmp88 + _tmp187 * _tmp298 - _tmp283 * _tmp44 * _tmp88);
  const Scalar _tmp306 = -_tmp287 * _tmp88 + _tmp298 * _tmp83;
  const Scalar _tmp307 = -_tmp126 * _tmp287 - _tmp145 * _tmp295 + _tmp146 * _tmp291 +
                         _tmp146 * _tmp292 - _tmp146 * _tmp294 + _tmp296 * _tmp83;
  const Scalar _tmp308 = _tmp208 * _tmp288;
  const Scalar _tmp309 = _tmp141 * _tmp293;
  const Scalar _tmp310 =
      -_tmp131 * _tmp287 + _tmp138 * _tmp308 - _tmp138 * _tmp309 + _tmp207 * _tmp285;
  const Scalar _tmp311 = -_tmp109 * _tmp287 - _tmp145 * _tmp303 - _tmp146 * _tmp300 +
                         _tmp146 * _tmp301 + _tmp146 * _tmp302 + _tmp299 * _tmp83;
  const Scalar _tmp312 = _tmp164 * (_tmp137 * (_tmp144 * (_tmp308 - _tmp309 - _tmp310 * _tmp94) +
                                               _tmp148 * (_tmp141 * _tmp291 + _tmp141 * _tmp292 -
                                                          _tmp141 * _tmp294 - _tmp307 * _tmp94) +
                                               _tmp154 * (-_tmp141 * _tmp300 + _tmp141 * _tmp301 +
                                                          _tmp141 * _tmp302 - _tmp311 * _tmp94) -
                                               _tmp203 * _tmp306) -
                                    _tmp201 * _tmp305);
  const Scalar _tmp313 = _tmp298 * _tmp90;
  const Scalar _tmp314 = _tmp144 * _tmp285;
  const Scalar _tmp315 = -_tmp221 * _tmp283 + _tmp222 * _tmp297 + _tmp223 * _tmp304 -
                         _tmp240 * _tmp283 + _tmp242 * _tmp313 + _tmp242 * _tmp314;
  const Scalar _tmp316 =
      _tmp247 *
      (_tmp225 * (_tmp226 * _tmp310 + _tmp228 * _tmp306 + _tmp229 * _tmp307 + _tmp230 * _tmp311) -
       _tmp248 * _tmp315);
  const Scalar _tmp317 = Scalar(9.6622558468725703) * _tmp315;
  const Scalar _tmp318 = -_tmp117 * _tmp299 + _tmp117 * _tmp303 * _tmp85 - _tmp129 * _tmp296 +
                         _tmp268 * _tmp295 - _tmp313 - _tmp314;
  const Scalar _tmp319 = Scalar(9.6622558468725703) * _tmp318;
  const Scalar _tmp320 =
      _tmp267 *
      (_tmp259 * (-_tmp134 * _tmp308 + _tmp134 * _tmp309 + _tmp260 * _tmp300 - _tmp260 * _tmp301 -
                  _tmp260 * _tmp302 - _tmp261 * _tmp291 - _tmp261 * _tmp292 + _tmp261 * _tmp294) -
       _tmp270 * _tmp318);

  // Output terms (1)
  Eigen::Matrix<Scalar, 4, 3> _res;

  _res(0, 0) = 0;
  _res(1, 0) =
      _tmp200 * _tmp219 +
      _tmp213 *
          (-_tmp162 * _tmp211 -
           _tmp218 * (-_tmp200 * _tmp216 + _tmp214 * (-_tmp200 * _tmp212 - _tmp211 * _tmp213)));
  _res(2, 0) =
      _tmp234 *
          (-_tmp249 * _tmp252 -
           _tmp251 * (_tmp236 * (-_tmp232 * _tmp245 - _tmp234 * _tmp249) - _tmp244 * _tmp250)) +
      _tmp238 * _tmp245;
  _res(3, 0) =
      -_tmp254 * _tmp257 +
      _tmp273 * (-_tmp265 * _tmp271 -
                 _tmp277 * (-_tmp269 * _tmp278 + _tmp275 * (-_tmp254 * _tmp279 - _tmp263 * _tmp280 -
                                                            _tmp271 * _tmp273))) +
      _tmp280 * _tmp281;
  _res(0, 1) = 0;
  _res(1, 1) =
      _tmp213 *
          (-_tmp162 * _tmp312 -
           _tmp218 * (_tmp214 * (-_tmp212 * _tmp305 - _tmp213 * _tmp312) - _tmp216 * _tmp305)) +
      _tmp219 * _tmp305;
  _res(2, 1) =
      _tmp234 *
          (-_tmp251 * (_tmp236 * (-_tmp232 * _tmp317 - _tmp234 * _tmp316) - _tmp250 * _tmp315) -
           _tmp252 * _tmp316) +
      _tmp238 * _tmp317;
  _res(3, 1) =
      -_tmp255 * _tmp257 +
      _tmp273 * (-_tmp265 * _tmp320 -
                 _tmp277 * (_tmp275 * (-_tmp255 * _tmp279 - _tmp263 * _tmp319 - _tmp273 * _tmp320) -
                            _tmp278 * _tmp318)) +
      _tmp281 * _tmp319;
  _res(0, 2) = 0;
  _res(1, 2) = 0;
  _res(2, 2) = 0;
  _res(3, 2) = -_tmp253 * _tmp257;

  return _res;
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym

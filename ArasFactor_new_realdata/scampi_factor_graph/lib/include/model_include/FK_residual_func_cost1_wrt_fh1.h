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
 * Symbolic function: FK_residual_func_cost1_wrt_fh1
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
 *     res: Matrix41
 */
template <typename Scalar>
Eigen::Matrix<Scalar, 4, 1> FkResidualFuncCost1WrtFh1(
    const Scalar fh1, const Scalar fv1, const sym::Rot3<Scalar>& DeltaRot,
    const sym::Pose3<Scalar>& TransformationMatrix, const Eigen::Matrix<Scalar, 4, 1>& encoder,
    const Eigen::Matrix<Scalar, 3, 1>& p_a, const Eigen::Matrix<Scalar, 3, 1>& p_b,
    const Eigen::Matrix<Scalar, 3, 1>& p_c, const Eigen::Matrix<Scalar, 3, 1>& p_d,
    const Scalar epsilon) {
  // Total ops: 639

  // Unused inputs
  (void)encoder;
  (void)epsilon;

  // Input arrays
  const Eigen::Matrix<Scalar, 4, 1>& _DeltaRot = DeltaRot.Data();
  const Eigen::Matrix<Scalar, 7, 1>& _TransformationMatrix = TransformationMatrix.Data();

  // Intermediate terms (218)
  const Scalar _tmp0 = Scalar(1.0) / (fh1);
  const Scalar _tmp1 =
      _DeltaRot[0] * _TransformationMatrix[3] - _DeltaRot[1] * _TransformationMatrix[2] +
      _DeltaRot[2] * _TransformationMatrix[1] + _DeltaRot[3] * _TransformationMatrix[0];
  const Scalar _tmp2 = -2 * std::pow(_tmp1, Scalar(2));
  const Scalar _tmp3 =
      -_DeltaRot[0] * _TransformationMatrix[1] + _DeltaRot[1] * _TransformationMatrix[0] +
      _DeltaRot[2] * _TransformationMatrix[3] + _DeltaRot[3] * _TransformationMatrix[2];
  const Scalar _tmp4 = -2 * std::pow(_tmp3, Scalar(2));
  const Scalar _tmp5 = Scalar(0.20999999999999999) * _tmp2 + Scalar(0.20999999999999999) * _tmp4 +
                       Scalar(0.20999999999999999);
  const Scalar _tmp6 = -_tmp5;
  const Scalar _tmp7 =
      _DeltaRot[0] * _TransformationMatrix[2] + _DeltaRot[1] * _TransformationMatrix[3] -
      _DeltaRot[2] * _TransformationMatrix[0] + _DeltaRot[3] * _TransformationMatrix[1];
  const Scalar _tmp8 = 2 * _tmp1;
  const Scalar _tmp9 = _tmp7 * _tmp8;
  const Scalar _tmp10 =
      -2 * _DeltaRot[0] * _TransformationMatrix[0] - 2 * _DeltaRot[1] * _TransformationMatrix[1] -
      2 * _DeltaRot[2] * _TransformationMatrix[2] + 2 * _DeltaRot[3] * _TransformationMatrix[3];
  const Scalar _tmp11 = _tmp10 * _tmp3;
  const Scalar _tmp12 = Scalar(0.20999999999999999) * _tmp11 + Scalar(0.20999999999999999) * _tmp9;
  const Scalar _tmp13 = 2 * _tmp3 * _tmp7;
  const Scalar _tmp14 = _tmp1 * _tmp10;
  const Scalar _tmp15 = _tmp13 - _tmp14;
  const Scalar _tmp16 = -Scalar(0.010999999999999999) * _tmp15;
  const Scalar _tmp17 = -_tmp12 + _tmp16;
  const Scalar _tmp18 = _tmp17 + _tmp6;
  const Scalar _tmp19 = _TransformationMatrix[5] + _tmp18;
  const Scalar _tmp20 = 1 - 2 * std::pow(_tmp7, Scalar(2));
  const Scalar _tmp21 = Scalar(0.20999999999999999) * _tmp20 + Scalar(0.20999999999999999) * _tmp4;
  const Scalar _tmp22 = -_tmp21;
  const Scalar _tmp23 = _tmp3 * _tmp8;
  const Scalar _tmp24 = _tmp10 * _tmp7;
  const Scalar _tmp25 = _tmp23 + _tmp24;
  const Scalar _tmp26 = -Scalar(0.010999999999999999) * _tmp25;
  const Scalar _tmp27 = -Scalar(0.20999999999999999) * _tmp11 + Scalar(0.20999999999999999) * _tmp9;
  const Scalar _tmp28 = _tmp26 - _tmp27;
  const Scalar _tmp29 = _tmp22 + _tmp28;
  const Scalar _tmp30 = _TransformationMatrix[4] + _tmp29;
  const Scalar _tmp31 = _tmp0 * fv1;
  const Scalar _tmp32 = std::asinh(_tmp31);
  const Scalar _tmp33 = Scalar(9.6622558468725703) * _tmp32;
  const Scalar _tmp34 =
      -Scalar(0.1034955) * _tmp33 * fh1 -
      Scalar(0.1034955) * std::sqrt(Scalar(std::pow(Scalar(-_tmp19 + p_a(1, 0)), Scalar(2)) +
                                           std::pow(Scalar(-_tmp30 + p_a(0, 0)), Scalar(2))));
  const Scalar _tmp35 = _tmp0 * _tmp34;
  const Scalar _tmp36 = Scalar(0.1034955) * p_a(2, 0);
  const Scalar _tmp37 = std::pow(fh1, Scalar(-2));
  const Scalar _tmp38 =
      std::pow(Scalar(_tmp37 * std::pow(fv1, Scalar(2)) + 1), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp39 = Scalar(1.0) * _tmp32;
  const Scalar _tmp40 = Scalar(43.164000000000001) - fv1;
  const Scalar _tmp41 = _tmp26 + _tmp27;
  const Scalar _tmp42 = _tmp21 + _tmp41;
  const Scalar _tmp43 = _TransformationMatrix[4] + _tmp42;
  const Scalar _tmp44 = _tmp43 - p_c(0, 0);
  const Scalar _tmp45 = _tmp12 + _tmp16;
  const Scalar _tmp46 = _tmp45 + _tmp5;
  const Scalar _tmp47 = _TransformationMatrix[5] + _tmp46;
  const Scalar _tmp48 = _tmp47 - p_c(1, 0);
  const Scalar _tmp49 = std::pow(Scalar(std::pow(_tmp44, Scalar(2)) + std::pow(_tmp48, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp50 = _tmp44 * _tmp49;
  const Scalar _tmp51 = _tmp45 + _tmp6;
  const Scalar _tmp52 = _TransformationMatrix[5] + _tmp51;
  const Scalar _tmp53 = _tmp52 - p_b(1, 0);
  const Scalar _tmp54 = _tmp21 + _tmp28;
  const Scalar _tmp55 = _TransformationMatrix[4] + _tmp54;
  const Scalar _tmp56 = _tmp55 - p_b(0, 0);
  const Scalar _tmp57 = Scalar(1.0) / (_tmp56);
  const Scalar _tmp58 = _tmp53 * _tmp57;
  const Scalar _tmp59 = _tmp48 * _tmp49;
  const Scalar _tmp60 = Scalar(1.0) / (_tmp50 * _tmp58 - _tmp59);
  const Scalar _tmp61 = Scalar(1.0) * _tmp54;
  const Scalar _tmp62 = Scalar(1.0) * _tmp51;
  const Scalar _tmp63 = -_tmp62;
  const Scalar _tmp64 = Scalar(1.0) / (_tmp46 + _tmp63);
  const Scalar _tmp65 = _tmp64 * (-_tmp42 + _tmp61);
  const Scalar _tmp66 = _tmp61 + _tmp62 * _tmp65;
  const Scalar _tmp67 = 0;
  const Scalar _tmp68 = Scalar(0.20999999999999999) * _tmp13 + Scalar(0.20999999999999999) * _tmp14;
  const Scalar _tmp69 =
      -Scalar(0.010999999999999999) * _tmp2 - Scalar(0.010999999999999999) * _tmp20;
  const Scalar _tmp70 = Scalar(0.20999999999999999) * _tmp23 - Scalar(0.20999999999999999) * _tmp24;
  const Scalar _tmp71 = _tmp69 - _tmp70;
  const Scalar _tmp72 = _tmp68 + _tmp71;
  const Scalar _tmp73 = _tmp17 + _tmp5;
  const Scalar _tmp74 = _TransformationMatrix[5] + _tmp73;
  const Scalar _tmp75 = _tmp74 - p_d(1, 0);
  const Scalar _tmp76 = _tmp22 + _tmp41;
  const Scalar _tmp77 = _TransformationMatrix[4] + _tmp76;
  const Scalar _tmp78 = _tmp77 - p_d(0, 0);
  const Scalar _tmp79 = std::pow(Scalar(std::pow(_tmp75, Scalar(2)) + std::pow(_tmp78, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp80 = _tmp78 * _tmp79;
  const Scalar _tmp81 = -_tmp68;
  const Scalar _tmp82 = _tmp69 + _tmp70;
  const Scalar _tmp83 = _tmp81 + _tmp82;
  const Scalar _tmp84 = _tmp68 + _tmp82;
  const Scalar _tmp85 = _tmp50 * _tmp83 - _tmp50 * _tmp84;
  const Scalar _tmp86 = _tmp75 * _tmp79;
  const Scalar _tmp87 = _tmp58 * _tmp80 - _tmp86;
  const Scalar _tmp88 = _tmp60 * _tmp87;
  const Scalar _tmp89 = _tmp58 * _tmp83;
  const Scalar _tmp90 = _tmp60 * (-_tmp50 * _tmp89 + _tmp59 * _tmp84);
  const Scalar _tmp91 = _tmp72 * _tmp86 - _tmp80 * _tmp89 - _tmp87 * _tmp90;
  const Scalar _tmp92 = -_tmp65 * _tmp91 - _tmp72 * _tmp80 + _tmp80 * _tmp83 - _tmp85 * _tmp88;
  const Scalar _tmp93 = Scalar(1.0) / (_tmp92);
  const Scalar _tmp94 = _tmp87 * _tmp93;
  const Scalar _tmp95 = _tmp60 * _tmp67 * _tmp94;
  const Scalar _tmp96 = _tmp67 * _tmp93;
  const Scalar _tmp97 =
      std::sqrt(Scalar(std::pow(_tmp53, Scalar(2)) + std::pow(_tmp56, Scalar(2))));
  const Scalar _tmp98 = _tmp57 * _tmp97;
  const Scalar _tmp99 = Scalar(1.0) * _tmp60;
  const Scalar _tmp100 = Scalar(1.0) * _tmp90;
  const Scalar _tmp101 = _tmp100 * _tmp65 - _tmp85 * _tmp99;
  const Scalar _tmp102 = Scalar(1.0) / (_tmp97);
  const Scalar _tmp103 = _tmp98 * (-_tmp102 * _tmp51 * _tmp56 + _tmp102 * _tmp53 * _tmp54);
  const Scalar _tmp104 = _tmp103 * _tmp50 - _tmp42 * _tmp59 + _tmp46 * _tmp50;
  const Scalar _tmp105 = _tmp103 * _tmp80 - _tmp104 * _tmp88 + _tmp73 * _tmp80 - _tmp76 * _tmp86;
  const Scalar _tmp106 = _tmp105 * _tmp93;
  const Scalar _tmp107 = Scalar(1.0) / (_tmp105);
  const Scalar _tmp108 = _tmp107 * _tmp92;
  const Scalar _tmp109 = _tmp108 * (-_tmp101 * _tmp106 - _tmp104 * _tmp99);
  const Scalar _tmp110 = _tmp101 + _tmp109;
  const Scalar _tmp111 = _tmp60 * (-_tmp110 * _tmp94 + Scalar(1.0));
  const Scalar _tmp112 = _tmp110 * _tmp93;
  const Scalar _tmp113 = _tmp30 - p_a(0, 0);
  const Scalar _tmp114 = _tmp19 - p_a(1, 0);
  const Scalar _tmp115 =
      std::pow(Scalar(std::pow(_tmp113, Scalar(2)) + std::pow(_tmp114, Scalar(2))),
               Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp116 = _tmp114 * _tmp115;
  const Scalar _tmp117 = _tmp116 * _tmp98 * (_tmp111 * _tmp50 + _tmp112 * _tmp80);
  const Scalar _tmp118 = _tmp58 * _tmp60;
  const Scalar _tmp119 = _tmp58 * _tmp90 + _tmp89;
  const Scalar _tmp120 = _tmp118 * _tmp85 - _tmp119 * _tmp65 - _tmp83;
  const Scalar _tmp121 = _tmp108 * (-_tmp103 + _tmp104 * _tmp118 - _tmp106 * _tmp120);
  const Scalar _tmp122 = _tmp120 + _tmp121;
  const Scalar _tmp123 = _tmp60 * (-_tmp122 * _tmp94 - _tmp58);
  const Scalar _tmp124 = _tmp122 * _tmp93;
  const Scalar _tmp125 = _tmp113 * _tmp115;
  const Scalar _tmp126 = _tmp125 * _tmp98 * (_tmp123 * _tmp50 + _tmp124 * _tmp80 + Scalar(1.0));
  const Scalar _tmp127 = _tmp116 * _tmp29 - _tmp125 * _tmp18;
  const Scalar _tmp128 = Scalar(1.0) * _tmp107;
  const Scalar _tmp129 = _tmp127 * _tmp98 * (-_tmp128 * _tmp50 * _tmp88 + _tmp128 * _tmp80);
  const Scalar _tmp130 = -_tmp117 * fh1 - _tmp126 * fh1 - _tmp129 * fh1 -
                         _tmp40 * _tmp98 * (-_tmp50 * _tmp95 + _tmp80 * _tmp96);
  const Scalar _tmp131 = Scalar(1.0) / (_tmp130);
  const Scalar _tmp132 = Scalar(0.1034955) * _tmp131;
  const Scalar _tmp133 = _tmp63 + _tmp73;
  const Scalar _tmp134 = _tmp133 * _tmp65;
  const Scalar _tmp135 = Scalar(1.0) / (-_tmp134 + _tmp61 - _tmp76);
  const Scalar _tmp136 = Scalar(1.0) * _tmp135;
  const Scalar _tmp137 = _tmp108 * _tmp136;
  const Scalar _tmp138 = _tmp64 * (-_tmp128 * _tmp91 + _tmp133 * _tmp137);
  const Scalar _tmp139 = Scalar(1.0) * _tmp127;
  const Scalar _tmp140 = _tmp139 * (_tmp137 - Scalar(1.0) * _tmp138);
  const Scalar _tmp141 = _tmp91 * _tmp93;
  const Scalar _tmp142 = _tmp133 * _tmp135;
  const Scalar _tmp143 = _tmp119 + _tmp121 * _tmp142 - _tmp122 * _tmp141;
  const Scalar _tmp144 = Scalar(1.0) * _tmp64;
  const Scalar _tmp145 = Scalar(1.0) * _tmp125 * (_tmp121 * _tmp136 - _tmp143 * _tmp144);
  const Scalar _tmp146 = _tmp71 + _tmp81;
  const Scalar _tmp147 = _tmp146 * fh1;
  const Scalar _tmp148 = -_tmp116 * _tmp147 - Scalar(5.1796800000000003) * _tmp15 - _tmp18 * fv1;
  const Scalar _tmp149 = _tmp136 * _tmp65;
  const Scalar _tmp150 = _tmp134 * _tmp136 + Scalar(1.0);
  const Scalar _tmp151 = -Scalar(1.0) * _tmp144 * _tmp150 + Scalar(1.0) * _tmp149;
  const Scalar _tmp152 = -_tmp100 + _tmp109 * _tmp142 - _tmp110 * _tmp141;
  const Scalar _tmp153 = Scalar(1.0) * _tmp116 * (_tmp109 * _tmp136 - _tmp144 * _tmp152);
  const Scalar _tmp154 = _tmp135 * _tmp66;
  const Scalar _tmp155 = -_tmp133 * _tmp154 - _tmp141 * _tmp67 + _tmp63;
  const Scalar _tmp156 = _tmp125 * _tmp147 + Scalar(5.1796800000000003) * _tmp25 + _tmp29 * fv1;
  const Scalar _tmp157 = _tmp133 * _tmp64;
  const Scalar _tmp158 = _tmp136 * _tmp157;
  const Scalar _tmp159 = -Scalar(1.0) * _tmp136 + Scalar(1.0) * _tmp158;
  const Scalar _tmp160 =
      _tmp140 * fh1 + _tmp145 * fh1 + _tmp148 * _tmp151 + _tmp153 * fh1 + _tmp156 * _tmp159 +
      Scalar(1.0) * _tmp40 * (-_tmp136 * _tmp66 - _tmp144 * _tmp155 + Scalar(1.0));
  const Scalar _tmp161 = std::asinh(_tmp131 * _tmp160);
  const Scalar _tmp162 = Scalar(9.6622558468725703) * _tmp130;
  const Scalar _tmp163 =
      -_tmp161 * _tmp162 - std::sqrt(Scalar(std::pow(Scalar(-_tmp52 + p_b(1, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp55 + p_b(0, 0)), Scalar(2))));
  const Scalar _tmp164 = _tmp132 * _tmp163;
  const Scalar _tmp165 = Scalar(1.0) * _tmp161;
  const Scalar _tmp166 = -_tmp117 - _tmp126 - _tmp129;
  const Scalar _tmp167 = Scalar(9.6622558468725703) * _tmp166;
  const Scalar _tmp168 = std::pow(_tmp130, Scalar(-2));
  const Scalar _tmp169 = _tmp166 * _tmp168;
  const Scalar _tmp170 = _tmp125 * _tmp146;
  const Scalar _tmp171 = _tmp116 * _tmp146;
  const Scalar _tmp172 =
      (_tmp131 * (_tmp140 + _tmp145 - _tmp151 * _tmp171 + _tmp153 + _tmp159 * _tmp170) -
       _tmp160 * _tmp169) /
      std::sqrt(Scalar(std::pow(_tmp160, Scalar(2)) * _tmp168 + 1));
  const Scalar _tmp173 = Scalar(0.1034955) * _tmp169;
  const Scalar _tmp174 = _tmp107 * _tmp139;
  const Scalar _tmp175 = _tmp174 * fh1;
  const Scalar _tmp176 = _tmp123 * _tmp125;
  const Scalar _tmp177 = _tmp111 * _tmp116;
  const Scalar _tmp178 = -_tmp175 * _tmp88 + _tmp176 * fh1 + _tmp177 * fh1 - _tmp40 * _tmp95;
  const Scalar _tmp179 = Scalar(1.0) / (_tmp178);
  const Scalar _tmp180 = _tmp150 * _tmp64;
  const Scalar _tmp181 = _tmp127 * _tmp138;
  const Scalar _tmp182 = _tmp136 * _tmp156;
  const Scalar _tmp183 = _tmp125 * _tmp143 * _tmp64;
  const Scalar _tmp184 = _tmp116 * _tmp152 * _tmp64;
  const Scalar _tmp185 = _tmp148 * _tmp180 + _tmp155 * _tmp40 * _tmp64 - _tmp157 * _tmp182 +
                         _tmp181 * fh1 + _tmp183 * fh1 + _tmp184 * fh1;
  const Scalar _tmp186 = std::asinh(_tmp179 * _tmp185);
  const Scalar _tmp187 = Scalar(1.0) * _tmp186;
  const Scalar _tmp188 = std::pow(_tmp178, Scalar(-2));
  const Scalar _tmp189 = -_tmp174 * _tmp88 + _tmp176 + _tmp177;
  const Scalar _tmp190 = _tmp188 * _tmp189;
  const Scalar _tmp191 =
      (_tmp179 * (-_tmp158 * _tmp170 - _tmp171 * _tmp180 + _tmp181 + _tmp183 + _tmp184) -
       _tmp185 * _tmp190) /
      std::sqrt(Scalar(std::pow(_tmp185, Scalar(2)) * _tmp188 + 1));
  const Scalar _tmp192 = Scalar(0.1034955) * _tmp190;
  const Scalar _tmp193 = Scalar(9.6622558468725703) * _tmp178;
  const Scalar _tmp194 =
      -_tmp186 * _tmp193 - std::sqrt(Scalar(std::pow(Scalar(-_tmp43 + p_c(0, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp47 + p_c(1, 0)), Scalar(2))));
  const Scalar _tmp195 = Scalar(9.6622558468725703) * _tmp189;
  const Scalar _tmp196 = Scalar(0.1034955) * _tmp179;
  const Scalar _tmp197 = _tmp194 * _tmp196;
  const Scalar _tmp198 = _tmp112 * _tmp116;
  const Scalar _tmp199 = _tmp124 * _tmp125;
  const Scalar _tmp200 = _tmp175 + _tmp198 * fh1 + _tmp199 * fh1 + _tmp40 * _tmp96;
  const Scalar _tmp201 = std::pow(_tmp200, Scalar(-2));
  const Scalar _tmp202 = _tmp174 + _tmp198 + _tmp199;
  const Scalar _tmp203 = _tmp201 * _tmp202;
  const Scalar _tmp204 = Scalar(0.1034955) * _tmp203;
  const Scalar _tmp205 = Scalar(1.0) / (_tmp200);
  const Scalar _tmp206 = _tmp109 * _tmp116 * _tmp135;
  const Scalar _tmp207 = _tmp121 * _tmp125 * _tmp135;
  const Scalar _tmp208 = _tmp127 * _tmp137;
  const Scalar _tmp209 = -_tmp148 * _tmp149 + _tmp154 * _tmp40 + _tmp182 - _tmp206 * fh1 -
                         _tmp207 * fh1 - _tmp208 * fh1;
  const Scalar _tmp210 = std::asinh(_tmp205 * _tmp209);
  const Scalar _tmp211 = Scalar(9.6622558468725703) * _tmp202;
  const Scalar _tmp212 = Scalar(9.6622558468725703) * _tmp200;
  const Scalar _tmp213 = (-_tmp203 * _tmp209 + _tmp205 * (_tmp136 * _tmp170 + _tmp149 * _tmp171 -
                                                          _tmp206 - _tmp207 - _tmp208)) /
                         std::sqrt(Scalar(_tmp201 * std::pow(_tmp209, Scalar(2)) + 1));
  const Scalar _tmp214 = Scalar(0.1034955) * _tmp205;
  const Scalar _tmp215 =
      -_tmp210 * _tmp212 - std::sqrt(Scalar(std::pow(Scalar(-_tmp74 + p_d(1, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp77 + p_d(0, 0)), Scalar(2))));
  const Scalar _tmp216 = _tmp214 * _tmp215;
  const Scalar _tmp217 = Scalar(1.0) * _tmp210;

  // Output terms (1)
  Eigen::Matrix<Scalar, 4, 1> _res;

  _res(0, 0) =
      -Scalar(9.6622558468725703) * _tmp0 * _tmp36 -
      Scalar(9.6622558468725703) * fh1 *
          (-_tmp36 * _tmp37 - Scalar(1.0) * _tmp37 * _tmp38 * fv1 * std::sinh(_tmp39) -
           (Scalar(0.1034955) * _tmp0 * (Scalar(9.6622558468725703) * _tmp31 * _tmp38 - _tmp33) -
            _tmp34 * _tmp37) *
               std::sinh(_tmp35)) +
      Scalar(9.6622558468725703) * std::cosh(_tmp35) -
      Scalar(9.6622558468725703) * std::cosh(_tmp39);
  _res(1, 0) =
      -_tmp162 * (Scalar(1.0) * _tmp172 * std::sinh(_tmp165) - _tmp173 * p_b(2, 0) -
                  (_tmp132 * (-_tmp161 * _tmp167 - _tmp162 * _tmp172) - _tmp163 * _tmp173) *
                      std::sinh(_tmp164)) -
      _tmp167 * (_tmp132 * p_b(2, 0) - std::cosh(_tmp164) + std::cosh(_tmp165));
  _res(2, 0) =
      -_tmp193 * (Scalar(1.0) * _tmp191 * std::sinh(_tmp187) - _tmp192 * p_c(2, 0) -
                  (-_tmp192 * _tmp194 + _tmp196 * (-_tmp186 * _tmp195 - _tmp191 * _tmp193)) *
                      std::sinh(_tmp197)) -
      _tmp195 * (_tmp196 * p_c(2, 0) + std::cosh(_tmp187) - std::cosh(_tmp197));
  _res(3, 0) =
      -_tmp211 * (_tmp214 * p_d(2, 0) - std::cosh(_tmp216) + std::cosh(_tmp217)) -
      _tmp212 * (-_tmp204 * p_d(2, 0) + Scalar(1.0) * _tmp213 * std::sinh(_tmp217) -
                 (-_tmp204 * _tmp215 + _tmp214 * (-_tmp210 * _tmp211 - _tmp212 * _tmp213)) *
                     std::sinh(_tmp216));

  return _res;
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym

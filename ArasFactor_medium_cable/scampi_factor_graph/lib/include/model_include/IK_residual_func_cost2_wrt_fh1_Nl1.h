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
 * Symbolic function: IK_residual_func_cost2_wrt_fh1_Nl1
 *
 * Args:
 *     fh1: Scalar
 *     fv1: Scalar
 *     DeltaRot: Rot3
 *     position_vector: Matrix31
 *     Rot_init: Rot3
 *     epsilon: Scalar
 *
 * Outputs:
 *     res: Matrix41
 */
template <typename Scalar>
Eigen::Matrix<Scalar, 4, 1> IkResidualFuncCost2WrtFh1Nl1(
    const Scalar fh1, const Scalar fv1, const sym::Rot3<Scalar>& DeltaRot,
    const Eigen::Matrix<Scalar, 3, 1>& position_vector, const sym::Rot3<Scalar>& Rot_init,
    const Scalar epsilon) {
  // Total ops: 637

  // Unused inputs
  (void)epsilon;

  // Input arrays
  const Eigen::Matrix<Scalar, 4, 1>& _DeltaRot = DeltaRot.Data();
  const Eigen::Matrix<Scalar, 4, 1>& _Rot_init = Rot_init.Data();

  // Intermediate terms (210)
  const Scalar _tmp0 = Scalar(1.0) / (fh1);
  const Scalar _tmp1 = _DeltaRot[0] * _Rot_init[3] - _DeltaRot[1] * _Rot_init[2] +
                       _DeltaRot[2] * _Rot_init[1] + _DeltaRot[3] * _Rot_init[0];
  const Scalar _tmp2 = _DeltaRot[0] * _Rot_init[2] + _DeltaRot[1] * _Rot_init[3] -
                       _DeltaRot[2] * _Rot_init[0] + _DeltaRot[3] * _Rot_init[1];
  const Scalar _tmp3 = 2 * _tmp2;
  const Scalar _tmp4 = _tmp1 * _tmp3;
  const Scalar _tmp5 = -_DeltaRot[0] * _Rot_init[1] + _DeltaRot[1] * _Rot_init[0] +
                       _DeltaRot[2] * _Rot_init[3] + _DeltaRot[3] * _Rot_init[2];
  const Scalar _tmp6 = -2 * _DeltaRot[0] * _Rot_init[0] - 2 * _DeltaRot[1] * _Rot_init[1] -
                       2 * _DeltaRot[2] * _Rot_init[2] + 2 * _DeltaRot[3] * _Rot_init[3];
  const Scalar _tmp7 = _tmp5 * _tmp6;
  const Scalar _tmp8 = Scalar(0.20999999999999999) * _tmp4 + Scalar(0.20999999999999999) * _tmp7;
  const Scalar _tmp9 = -_tmp8;
  const Scalar _tmp10 = -2 * std::pow(_tmp1, Scalar(2));
  const Scalar _tmp11 = -2 * std::pow(_tmp5, Scalar(2));
  const Scalar _tmp12 = Scalar(0.20999999999999999) * _tmp10 +
                        Scalar(0.20999999999999999) * _tmp11 + Scalar(0.20999999999999999);
  const Scalar _tmp13 = _tmp3 * _tmp5;
  const Scalar _tmp14 = _tmp1 * _tmp6;
  const Scalar _tmp15 = _tmp13 - _tmp14;
  const Scalar _tmp16 = -Scalar(0.010999999999999999) * _tmp15;
  const Scalar _tmp17 = -_tmp12 + _tmp16;
  const Scalar _tmp18 = _tmp17 + _tmp9;
  const Scalar _tmp19 = _tmp18 + position_vector(1, 0);
  const Scalar _tmp20 = Scalar(0.20999999999999999) * _tmp4 - Scalar(0.20999999999999999) * _tmp7;
  const Scalar _tmp21 = -_tmp20;
  const Scalar _tmp22 = 2 * _tmp1 * _tmp5;
  const Scalar _tmp23 = _tmp2 * _tmp6;
  const Scalar _tmp24 = _tmp22 + _tmp23;
  const Scalar _tmp25 = -Scalar(0.010999999999999999) * _tmp24;
  const Scalar _tmp26 = 1 - 2 * std::pow(_tmp2, Scalar(2));
  const Scalar _tmp27 = Scalar(0.20999999999999999) * _tmp11 + Scalar(0.20999999999999999) * _tmp26;
  const Scalar _tmp28 = _tmp25 - _tmp27;
  const Scalar _tmp29 = _tmp21 + _tmp28;
  const Scalar _tmp30 = _tmp29 + position_vector(0, 0);
  const Scalar _tmp31 = _tmp0 * fv1;
  const Scalar _tmp32 = std::asinh(_tmp31);
  const Scalar _tmp33 = Scalar(9.6622558468725703) * _tmp32;
  const Scalar _tmp34 =
      -Scalar(0.1034955) * _tmp33 * fh1 -
      Scalar(0.86104699605560286) *
          std::sqrt(
              Scalar(std::pow(Scalar(-Scalar(0.12019727201198747) * _tmp19 - 1), Scalar(2)) +
                     Scalar(0.057067943527034905) *
                         std::pow(Scalar(-Scalar(0.50315118477277876) * _tmp30 - 1), Scalar(2))));
  const Scalar _tmp35 = _tmp0 * _tmp34;
  const Scalar _tmp36 = std::pow(fh1, Scalar(-2));
  const Scalar _tmp37 =
      std::pow(Scalar(_tmp36 * std::pow(fv1, Scalar(2)) + 1), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp38 = Scalar(1.0) * _tmp32;
  const Scalar _tmp39 = Scalar(43.164000000000001) - fv1;
  const Scalar _tmp40 = _tmp25 + _tmp27;
  const Scalar _tmp41 = _tmp21 + _tmp40;
  const Scalar _tmp42 = Scalar(1.0) * _tmp41;
  const Scalar _tmp43 = _tmp17 + _tmp8;
  const Scalar _tmp44 = Scalar(1.0) * _tmp43;
  const Scalar _tmp45 = _tmp20 + _tmp28;
  const Scalar _tmp46 = -_tmp44;
  const Scalar _tmp47 = _tmp12 + _tmp16;
  const Scalar _tmp48 = _tmp47 + _tmp9;
  const Scalar _tmp49 = Scalar(1.0) / (_tmp46 + _tmp48);
  const Scalar _tmp50 = _tmp49 * (_tmp42 - _tmp45);
  const Scalar _tmp51 = _tmp42 + _tmp44 * _tmp50;
  const Scalar _tmp52 = _tmp47 + _tmp8;
  const Scalar _tmp53 = _tmp46 + _tmp52;
  const Scalar _tmp54 = _tmp50 * _tmp53;
  const Scalar _tmp55 = _tmp20 + _tmp40;
  const Scalar _tmp56 = Scalar(1.0) / (_tmp42 - _tmp54 - _tmp55);
  const Scalar _tmp57 = Scalar(1.0) * _tmp56;
  const Scalar _tmp58 = _tmp52 + position_vector(1, 0);
  const Scalar _tmp59 = _tmp58 + Scalar(-4.7744369927459998);
  const Scalar _tmp60 = _tmp55 + position_vector(0, 0);
  const Scalar _tmp61 = _tmp60 + Scalar(-2.7171519410699099);
  const Scalar _tmp62 = std::pow(Scalar(std::pow(_tmp59, Scalar(2)) + std::pow(_tmp61, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp63 = _tmp61 * _tmp62;
  const Scalar _tmp64 = Scalar(0.20999999999999999) * _tmp13 + Scalar(0.20999999999999999) * _tmp14;
  const Scalar _tmp65 = -_tmp64;
  const Scalar _tmp66 =
      -Scalar(0.010999999999999999) * _tmp10 - Scalar(0.010999999999999999) * _tmp26;
  const Scalar _tmp67 = Scalar(0.20999999999999999) * _tmp22 - Scalar(0.20999999999999999) * _tmp23;
  const Scalar _tmp68 = _tmp66 + _tmp67;
  const Scalar _tmp69 = _tmp65 + _tmp68;
  const Scalar _tmp70 = _tmp41 + position_vector(0, 0);
  const Scalar _tmp71 = _tmp70 + Scalar(-2.5193355532036801);
  const Scalar _tmp72 = Scalar(1.0) / (_tmp71);
  const Scalar _tmp73 = _tmp43 + position_vector(1, 0);
  const Scalar _tmp74 = _tmp73 + Scalar(8.3885017487099702);
  const Scalar _tmp75 = _tmp72 * _tmp74;
  const Scalar _tmp76 = _tmp69 * _tmp75;
  const Scalar _tmp77 = _tmp64 + _tmp68;
  const Scalar _tmp78 = _tmp59 * _tmp62;
  const Scalar _tmp79 = _tmp45 + position_vector(0, 0);
  const Scalar _tmp80 = _tmp79 + Scalar(1.7965602546229);
  const Scalar _tmp81 = _tmp48 + position_vector(1, 0);
  const Scalar _tmp82 = _tmp81 + Scalar(-4.83288938413423);
  const Scalar _tmp83 = std::pow(Scalar(std::pow(_tmp80, Scalar(2)) + std::pow(_tmp82, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp84 = _tmp80 * _tmp83;
  const Scalar _tmp85 = _tmp69 * _tmp84;
  const Scalar _tmp86 = _tmp66 - _tmp67;
  const Scalar _tmp87 = _tmp64 + _tmp86;
  const Scalar _tmp88 = _tmp82 * _tmp83;
  const Scalar _tmp89 = -_tmp75 * _tmp85 + _tmp87 * _tmp88;
  const Scalar _tmp90 = _tmp63 * _tmp75 - _tmp78;
  const Scalar _tmp91 = Scalar(1.0) / (_tmp75 * _tmp84 - _tmp88);
  const Scalar _tmp92 = _tmp90 * _tmp91;
  const Scalar _tmp93 = -_tmp63 * _tmp76 + _tmp77 * _tmp78 - _tmp89 * _tmp92;
  const Scalar _tmp94 = -_tmp84 * _tmp87 + _tmp85;
  const Scalar _tmp95 = -_tmp50 * _tmp93 + _tmp63 * _tmp69 - _tmp63 * _tmp77 - _tmp92 * _tmp94;
  const Scalar _tmp96 = Scalar(1.0) / (_tmp95);
  const Scalar _tmp97 = 0;
  const Scalar _tmp98 = _tmp51 * _tmp56;
  const Scalar _tmp99 = _tmp46 - _tmp53 * _tmp98 - _tmp93 * _tmp97;
  const Scalar _tmp100 = Scalar(1.0) * _tmp49;
  const Scalar _tmp101 = _tmp30 + Scalar(1.9874742031097401);
  const Scalar _tmp102 = _tmp19 + Scalar(8.3196563720703107);
  const Scalar _tmp103 =
      std::pow(Scalar(std::pow(_tmp101, Scalar(2)) + std::pow(_tmp102, Scalar(2))),
               Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp104 = _tmp102 * _tmp103;
  const Scalar _tmp105 = _tmp65 + _tmp86;
  const Scalar _tmp106 = _tmp105 * fh1;
  const Scalar _tmp107 = -_tmp104 * _tmp106 - Scalar(5.1796800000000003) * _tmp15 - _tmp18 * fv1;
  const Scalar _tmp108 = _tmp54 * _tmp57 + Scalar(1.0);
  const Scalar _tmp109 = _tmp50 * _tmp57;
  const Scalar _tmp110 = -Scalar(1.0) * _tmp100 * _tmp108 + Scalar(1.0) * _tmp109;
  const Scalar _tmp111 =
      std::sqrt(Scalar(std::pow(_tmp71, Scalar(2)) + std::pow(_tmp74, Scalar(2))));
  const Scalar _tmp112 = Scalar(1.0) / (_tmp111);
  const Scalar _tmp113 = _tmp111 * _tmp72;
  const Scalar _tmp114 = _tmp113 * (_tmp112 * _tmp41 * _tmp74 - _tmp112 * _tmp43 * _tmp71);
  const Scalar _tmp115 = _tmp114 * _tmp84 - _tmp45 * _tmp88 + _tmp48 * _tmp84;
  const Scalar _tmp116 = _tmp114 * _tmp63 - _tmp115 * _tmp92 + _tmp52 * _tmp63 - _tmp55 * _tmp78;
  const Scalar _tmp117 = Scalar(1.0) / (_tmp116);
  const Scalar _tmp118 = Scalar(1.0) * _tmp117;
  const Scalar _tmp119 = _tmp117 * _tmp95;
  const Scalar _tmp120 = _tmp119 * _tmp57;
  const Scalar _tmp121 = -_tmp118 * _tmp93 + _tmp120 * _tmp53;
  const Scalar _tmp122 = _tmp101 * _tmp103;
  const Scalar _tmp123 = _tmp104 * _tmp29 - _tmp122 * _tmp18;
  const Scalar _tmp124 = Scalar(1.0) * _tmp123 * (-_tmp100 * _tmp121 + _tmp120);
  const Scalar _tmp125 = _tmp106 * _tmp122 + Scalar(5.1796800000000003) * _tmp24 + _tmp29 * fv1;
  const Scalar _tmp126 = _tmp49 * _tmp53;
  const Scalar _tmp127 = _tmp126 * _tmp57;
  const Scalar _tmp128 = Scalar(1.0) * _tmp127 - Scalar(1.0) * _tmp57;
  const Scalar _tmp129 = Scalar(1.0) * _tmp91;
  const Scalar _tmp130 = _tmp129 * _tmp89;
  const Scalar _tmp131 = -_tmp129 * _tmp94 + _tmp130 * _tmp50;
  const Scalar _tmp132 = _tmp116 * _tmp96;
  const Scalar _tmp133 = _tmp119 * (-_tmp115 * _tmp129 - _tmp131 * _tmp132);
  const Scalar _tmp134 = _tmp96 * (_tmp131 + _tmp133);
  const Scalar _tmp135 = _tmp53 * _tmp56;
  const Scalar _tmp136 = -_tmp130 + _tmp133 * _tmp135 - _tmp134 * _tmp93;
  const Scalar _tmp137 = Scalar(1.0) * _tmp104 * (-_tmp100 * _tmp136 + _tmp133 * _tmp57);
  const Scalar _tmp138 = _tmp75 * _tmp91;
  const Scalar _tmp139 = _tmp138 * _tmp89 + _tmp76;
  const Scalar _tmp140 = _tmp138 * _tmp94 - _tmp139 * _tmp50 - _tmp69;
  const Scalar _tmp141 = _tmp119 * (-_tmp114 + _tmp115 * _tmp138 - _tmp132 * _tmp140);
  const Scalar _tmp142 = _tmp96 * (_tmp140 + _tmp141);
  const Scalar _tmp143 = _tmp135 * _tmp141 + _tmp139 - _tmp142 * _tmp93;
  const Scalar _tmp144 = Scalar(1.0) * _tmp122 * (-_tmp100 * _tmp143 + _tmp141 * _tmp57);
  const Scalar _tmp145 = _tmp107 * _tmp110 + _tmp124 * fh1 + _tmp125 * _tmp128 + _tmp137 * fh1 +
                         _tmp144 * fh1 +
                         Scalar(1.0) * _tmp39 * (-_tmp100 * _tmp99 - _tmp51 * _tmp57 + Scalar(1.0));
  const Scalar _tmp146 = _tmp84 * _tmp92;
  const Scalar _tmp147 = _tmp113 * _tmp123 * (-_tmp118 * _tmp146 + _tmp118 * _tmp63);
  const Scalar _tmp148 = _tmp91 * (-_tmp134 * _tmp90 + Scalar(1.0));
  const Scalar _tmp149 = _tmp104 * _tmp113 * (_tmp134 * _tmp63 + _tmp148 * _tmp84);
  const Scalar _tmp150 = _tmp91 * (-_tmp142 * _tmp90 - _tmp75);
  const Scalar _tmp151 = _tmp113 * _tmp122 * (_tmp142 * _tmp63 + _tmp150 * _tmp84 + Scalar(1.0));
  const Scalar _tmp152 = -_tmp113 * _tmp39 * (-_tmp146 * _tmp97 + _tmp63 * _tmp97) - _tmp147 * fh1 -
                         _tmp149 * fh1 - _tmp151 * fh1;
  const Scalar _tmp153 = Scalar(1.0) / (_tmp152);
  const Scalar _tmp154 = std::asinh(_tmp145 * _tmp153);
  const Scalar _tmp155 = -_tmp147 - _tmp149 - _tmp151;
  const Scalar _tmp156 = Scalar(9.6622558468725703) * _tmp155;
  const Scalar _tmp157 = Scalar(9.6622558468725703) * _tmp152;
  const Scalar _tmp158 = _tmp105 * _tmp122;
  const Scalar _tmp159 = _tmp104 * _tmp105;
  const Scalar _tmp160 = std::pow(_tmp152, Scalar(-2));
  const Scalar _tmp161 = _tmp155 * _tmp160;
  const Scalar _tmp162 = (-_tmp145 * _tmp161 + _tmp153 * (-_tmp110 * _tmp159 + _tmp124 +
                                                          _tmp128 * _tmp158 + _tmp137 + _tmp144)) /
                         std::sqrt(Scalar(std::pow(_tmp145, Scalar(2)) * _tmp160 + 1));
  const Scalar _tmp163 = Scalar(0.1034955) * _tmp153;
  const Scalar _tmp164 =
      -_tmp154 * _tmp157 -
      Scalar(8.3885017487099702) *
          std::sqrt(
              Scalar(Scalar(0.090199313518583735) *
                         std::pow(Scalar(1 - Scalar(0.39693005512043167) * _tmp70), Scalar(2)) +
                     std::pow(Scalar(-Scalar(0.11921079949155229) * _tmp73 - 1), Scalar(2))));
  const Scalar _tmp165 = _tmp163 * _tmp164;
  const Scalar _tmp166 = Scalar(1.0) * _tmp154;
  const Scalar _tmp167 = _tmp122 * _tmp150;
  const Scalar _tmp168 = _tmp104 * _tmp148;
  const Scalar _tmp169 = _tmp39 * _tmp97;
  const Scalar _tmp170 = _tmp118 * _tmp123;
  const Scalar _tmp171 = _tmp170 * fh1;
  const Scalar _tmp172 = _tmp167 * fh1 + _tmp168 * fh1 - _tmp169 * _tmp92 - _tmp171 * _tmp92;
  const Scalar _tmp173 = Scalar(1.0) / (_tmp172);
  const Scalar _tmp174 = _tmp121 * _tmp123 * _tmp49;
  const Scalar _tmp175 = _tmp122 * _tmp143 * _tmp49;
  const Scalar _tmp176 = _tmp108 * _tmp49;
  const Scalar _tmp177 = _tmp125 * _tmp57;
  const Scalar _tmp178 = _tmp104 * _tmp136 * _tmp49;
  const Scalar _tmp179 = _tmp107 * _tmp176 - _tmp126 * _tmp177 + _tmp174 * fh1 + _tmp175 * fh1 +
                         _tmp178 * fh1 + _tmp39 * _tmp49 * _tmp99;
  const Scalar _tmp180 = std::asinh(_tmp173 * _tmp179);
  const Scalar _tmp181 = Scalar(9.6622558468725703) * _tmp172;
  const Scalar _tmp182 =
      -_tmp180 * _tmp181 -
      Scalar(4.83288938413423) *
          std::sqrt(
              Scalar(std::pow(Scalar(1 - Scalar(0.20691555724053495) * _tmp81), Scalar(2)) +
                     Scalar(0.13818785160942856) *
                         std::pow(Scalar(-Scalar(0.55661923802822921) * _tmp79 - 1), Scalar(2))));
  const Scalar _tmp183 = Scalar(0.1034955) * _tmp173;
  const Scalar _tmp184 = _tmp182 * _tmp183;
  const Scalar _tmp185 = Scalar(1.0) * _tmp180;
  const Scalar _tmp186 = _tmp167 + _tmp168 - _tmp170 * _tmp92;
  const Scalar _tmp187 = Scalar(9.6622558468725703) * _tmp186;
  const Scalar _tmp188 = std::pow(_tmp172, Scalar(-2));
  const Scalar _tmp189 = _tmp186 * _tmp188;
  const Scalar _tmp190 =
      (_tmp173 * (-_tmp127 * _tmp158 - _tmp159 * _tmp176 + _tmp174 + _tmp175 + _tmp178) -
       _tmp179 * _tmp189) /
      std::sqrt(Scalar(std::pow(_tmp179, Scalar(2)) * _tmp188 + 1));
  const Scalar _tmp191 = _tmp104 * _tmp134;
  const Scalar _tmp192 = _tmp122 * _tmp142;
  const Scalar _tmp193 = _tmp169 + _tmp171 + _tmp191 * fh1 + _tmp192 * fh1;
  const Scalar _tmp194 = _tmp122 * _tmp141 * _tmp56;
  const Scalar _tmp195 = _tmp104 * _tmp133 * _tmp56;
  const Scalar _tmp196 = _tmp120 * _tmp123;
  const Scalar _tmp197 = -_tmp107 * _tmp109 + _tmp177 - _tmp194 * fh1 - _tmp195 * fh1 -
                         _tmp196 * fh1 + _tmp39 * _tmp98;
  const Scalar _tmp198 = Scalar(1.0) / (_tmp193);
  const Scalar _tmp199 = std::asinh(_tmp197 * _tmp198);
  const Scalar _tmp200 = Scalar(9.6622558468725703) * _tmp199;
  const Scalar _tmp201 =
      -_tmp193 * _tmp200 -
      Scalar(4.7744369927459998) *
          std::sqrt(
              Scalar(std::pow(Scalar(1 - Scalar(0.2094487793051498) * _tmp58), Scalar(2)) +
                     Scalar(0.32387954179207445) *
                         std::pow(Scalar(1 - Scalar(0.36803241838814449) * _tmp60), Scalar(2))));
  const Scalar _tmp202 = Scalar(0.1034955) * _tmp198;
  const Scalar _tmp203 = _tmp201 * _tmp202;
  const Scalar _tmp204 = Scalar(1.0) * _tmp199;
  const Scalar _tmp205 = _tmp170 + _tmp191 + _tmp192;
  const Scalar _tmp206 = Scalar(9.6622558468725703) * _tmp193;
  const Scalar _tmp207 = std::pow(_tmp193, Scalar(-2));
  const Scalar _tmp208 = _tmp205 * _tmp207;
  const Scalar _tmp209 = (-_tmp197 * _tmp208 + _tmp198 * (_tmp109 * _tmp159 + _tmp158 * _tmp57 -
                                                          _tmp194 - _tmp195 - _tmp196)) /
                         std::sqrt(Scalar(std::pow(_tmp197, Scalar(2)) * _tmp207 + 1));

  // Output terms (1)
  Eigen::Matrix<Scalar, 4, 1> _res;

  _res(0, 0) =
      Scalar(9.6622558468725703) * fh1 *
          (Scalar(1.0) * _tmp36 * _tmp37 * fv1 * std::cosh(_tmp38) -
           (Scalar(0.1034955) * _tmp0 * (Scalar(9.6622558468725703) * _tmp31 * _tmp37 - _tmp33) -
            _tmp34 * _tmp36) *
               std::cosh(_tmp35)) -
      Scalar(9.6622558468725703) * std::sinh(_tmp35) -
      Scalar(9.6622558468725703) * std::sinh(_tmp38);
  _res(1, 0) = _tmp156 * (-std::sinh(_tmp165) - std::sinh(_tmp166)) +
               _tmp157 * (-Scalar(1.0) * _tmp162 * std::cosh(_tmp166) -
                          (-Scalar(0.1034955) * _tmp161 * _tmp164 +
                           _tmp163 * (-_tmp154 * _tmp156 - _tmp157 * _tmp162)) *
                              std::cosh(_tmp165));
  _res(2, 0) = _tmp181 * (-Scalar(1.0) * _tmp190 * std::cosh(_tmp185) -
                          (-Scalar(0.1034955) * _tmp182 * _tmp189 +
                           _tmp183 * (-_tmp180 * _tmp187 - _tmp181 * _tmp190)) *
                              std::cosh(_tmp184)) +
               _tmp187 * (-std::sinh(_tmp184) - std::sinh(_tmp185));
  _res(3, 0) = Scalar(9.6622558468725703) * _tmp205 * (-std::sinh(_tmp203) - std::sinh(_tmp204)) +
               _tmp206 * (-Scalar(1.0) * _tmp209 * std::cosh(_tmp204) -
                          (-Scalar(0.1034955) * _tmp201 * _tmp208 +
                           _tmp202 * (-_tmp200 * _tmp205 - _tmp206 * _tmp209)) *
                              std::cosh(_tmp203));

  return _res;
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym

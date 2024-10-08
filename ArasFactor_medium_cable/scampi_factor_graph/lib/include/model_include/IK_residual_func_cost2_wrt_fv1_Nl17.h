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
 * Symbolic function: IK_residual_func_cost2_wrt_fv1_Nl17
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
Eigen::Matrix<Scalar, 4, 1> IkResidualFuncCost2WrtFv1Nl17(
    const Scalar fh1, const Scalar fv1, const sym::Rot3<Scalar>& DeltaRot,
    const Eigen::Matrix<Scalar, 3, 1>& position_vector, const sym::Rot3<Scalar>& Rot_init,
    const Scalar epsilon) {
  // Total ops: 598

  // Unused inputs
  (void)epsilon;

  // Input arrays
  const Eigen::Matrix<Scalar, 4, 1>& _DeltaRot = DeltaRot.Data();
  const Eigen::Matrix<Scalar, 4, 1>& _Rot_init = Rot_init.Data();

  // Intermediate terms (191)
  const Scalar _tmp0 = Scalar(1.0) / (fh1);
  const Scalar _tmp1 = std::asinh(_tmp0 * fv1);
  const Scalar _tmp2 = Scalar(1.0) * _tmp0 /
                       std::sqrt(Scalar(1 + std::pow(fv1, Scalar(2)) / std::pow(fh1, Scalar(2))));
  const Scalar _tmp3 = _DeltaRot[0] * _Rot_init[3] - _DeltaRot[1] * _Rot_init[2] +
                       _DeltaRot[2] * _Rot_init[1] + _DeltaRot[3] * _Rot_init[0];
  const Scalar _tmp4 = _DeltaRot[0] * _Rot_init[2] + _DeltaRot[1] * _Rot_init[3] -
                       _DeltaRot[2] * _Rot_init[0] + _DeltaRot[3] * _Rot_init[1];
  const Scalar _tmp5 = 2 * _tmp4;
  const Scalar _tmp6 = _tmp3 * _tmp5;
  const Scalar _tmp7 = -_DeltaRot[0] * _Rot_init[1] + _DeltaRot[1] * _Rot_init[0] +
                       _DeltaRot[2] * _Rot_init[3] + _DeltaRot[3] * _Rot_init[2];
  const Scalar _tmp8 = -2 * _DeltaRot[0] * _Rot_init[0] - 2 * _DeltaRot[1] * _Rot_init[1] -
                       2 * _DeltaRot[2] * _Rot_init[2] + 2 * _DeltaRot[3] * _Rot_init[3];
  const Scalar _tmp9 = _tmp7 * _tmp8;
  const Scalar _tmp10 = Scalar(0.20999999999999999) * _tmp6 - Scalar(0.20999999999999999) * _tmp9;
  const Scalar _tmp11 = -2 * std::pow(_tmp4, Scalar(2));
  const Scalar _tmp12 = 1 - 2 * std::pow(_tmp7, Scalar(2));
  const Scalar _tmp13 = Scalar(0.20999999999999999) * _tmp11 + Scalar(0.20999999999999999) * _tmp12;
  const Scalar _tmp14 = 2 * _tmp3 * _tmp7;
  const Scalar _tmp15 = _tmp4 * _tmp8;
  const Scalar _tmp16 = _tmp14 + _tmp15;
  const Scalar _tmp17 = -Scalar(0.010999999999999999) * _tmp16;
  const Scalar _tmp18 = _tmp13 + _tmp17;
  const Scalar _tmp19 = _tmp10 + _tmp18;
  const Scalar _tmp20 = _tmp19 + position_vector(0, 0);
  const Scalar _tmp21 = Scalar(0.20999999999999999) * _tmp6 + Scalar(0.20999999999999999) * _tmp9;
  const Scalar _tmp22 = -2 * std::pow(_tmp3, Scalar(2));
  const Scalar _tmp23 = Scalar(0.20999999999999999) * _tmp12 + Scalar(0.20999999999999999) * _tmp22;
  const Scalar _tmp24 = _tmp5 * _tmp7;
  const Scalar _tmp25 = _tmp3 * _tmp8;
  const Scalar _tmp26 = _tmp24 - _tmp25;
  const Scalar _tmp27 = Scalar(0.010999999999999999) * _tmp26;
  const Scalar _tmp28 = -_tmp27;
  const Scalar _tmp29 = _tmp23 + _tmp28;
  const Scalar _tmp30 = _tmp21 + _tmp29;
  const Scalar _tmp31 = _tmp30 + position_vector(1, 0);
  const Scalar _tmp32 = Scalar(9.6622558468725703) * fh1;
  const Scalar _tmp33 = Scalar(43.164000000000001) - fv1;
  const Scalar _tmp34 = Scalar(0.20999999999999999) * _tmp24 + Scalar(0.20999999999999999) * _tmp25;
  const Scalar _tmp35 = -_tmp34;
  const Scalar _tmp36 = -Scalar(0.010999999999999999) * _tmp11 -
                        Scalar(0.010999999999999999) * _tmp22 + Scalar(-0.010999999999999999);
  const Scalar _tmp37 = Scalar(0.20999999999999999) * _tmp14 - Scalar(0.20999999999999999) * _tmp15;
  const Scalar _tmp38 = _tmp36 - _tmp37;
  const Scalar _tmp39 = _tmp35 + _tmp38;
  const Scalar _tmp40 = -_tmp10;
  const Scalar _tmp41 = -_tmp13 + _tmp17;
  const Scalar _tmp42 = _tmp40 + _tmp41;
  const Scalar _tmp43 = _tmp42 + position_vector(0, 0);
  const Scalar _tmp44 = _tmp43 + Scalar(1.9874742031097401);
  const Scalar _tmp45 = -_tmp21;
  const Scalar _tmp46 = -_tmp23;
  const Scalar _tmp47 = _tmp28 + _tmp46;
  const Scalar _tmp48 = _tmp45 + _tmp47;
  const Scalar _tmp49 = _tmp48 + position_vector(1, 0);
  const Scalar _tmp50 = _tmp49 + Scalar(8.3196563720703107);
  const Scalar _tmp51 = std::pow(Scalar(std::pow(_tmp44, Scalar(2)) + std::pow(_tmp50, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp52 = _tmp44 * _tmp51;
  const Scalar _tmp53 = _tmp50 * _tmp51;
  const Scalar _tmp54 = _tmp29 + _tmp45;
  const Scalar _tmp55 = _tmp54 + position_vector(1, 0);
  const Scalar _tmp56 = _tmp55 + Scalar(-4.83288938413423);
  const Scalar _tmp57 = _tmp10 + _tmp41;
  const Scalar _tmp58 = _tmp57 + position_vector(0, 0);
  const Scalar _tmp59 = _tmp58 + Scalar(1.7965602546229);
  const Scalar _tmp60 = Scalar(1.0) / (_tmp59);
  const Scalar _tmp61 = _tmp56 * _tmp60;
  const Scalar _tmp62 = _tmp52 * _tmp61 - _tmp53;
  const Scalar _tmp63 = _tmp21 + _tmp47;
  const Scalar _tmp64 = _tmp63 + position_vector(1, 0);
  const Scalar _tmp65 = _tmp64 + Scalar(8.3885017487099702);
  const Scalar _tmp66 = _tmp18 + _tmp40;
  const Scalar _tmp67 = _tmp66 + position_vector(0, 0);
  const Scalar _tmp68 = _tmp67 + Scalar(-2.5193355532036801);
  const Scalar _tmp69 = std::pow(Scalar(std::pow(_tmp65, Scalar(2)) + std::pow(_tmp68, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp70 = _tmp68 * _tmp69;
  const Scalar _tmp71 = _tmp65 * _tmp69;
  const Scalar _tmp72 = Scalar(1.0) / (_tmp61 * _tmp70 - _tmp71);
  const Scalar _tmp73 = _tmp36 + _tmp37;
  const Scalar _tmp74 = _tmp35 + _tmp73;
  const Scalar _tmp75 = _tmp34 + _tmp38;
  const Scalar _tmp76 = _tmp70 * _tmp75;
  const Scalar _tmp77 = _tmp72 * (-_tmp61 * _tmp76 + _tmp71 * _tmp74);
  const Scalar _tmp78 = _tmp61 * _tmp75;
  const Scalar _tmp79 = _tmp39 * _tmp53 - _tmp52 * _tmp78 - _tmp62 * _tmp77;
  const Scalar _tmp80 = Scalar(1.0) * _tmp54;
  const Scalar _tmp81 = -_tmp80;
  const Scalar _tmp82 = Scalar(1.0) / (_tmp63 + _tmp81);
  const Scalar _tmp83 = Scalar(1.0) * _tmp57;
  const Scalar _tmp84 = _tmp82 * (-_tmp66 + _tmp83);
  const Scalar _tmp85 = _tmp72 * (-_tmp70 * _tmp74 + _tmp76);
  const Scalar _tmp86 = -_tmp39 * _tmp52 + _tmp52 * _tmp75 - _tmp62 * _tmp85 - _tmp79 * _tmp84;
  const Scalar _tmp87 = Scalar(1.0) / (_tmp86);
  const Scalar _tmp88 = _tmp80 * _tmp84 + _tmp83;
  const Scalar _tmp89 = 0;
  const Scalar _tmp90 = _tmp70 * _tmp72;
  const Scalar _tmp91 = _tmp62 * _tmp90;
  const Scalar _tmp92 =
      std::sqrt(Scalar(std::pow(_tmp56, Scalar(2)) + std::pow(_tmp59, Scalar(2))));
  const Scalar _tmp93 = _tmp60 * _tmp92;
  const Scalar _tmp94 = _tmp93 * (_tmp52 * _tmp89 - _tmp89 * _tmp91);
  const Scalar _tmp95 = Scalar(1.0) / (_tmp92);
  const Scalar _tmp96 = _tmp93 * (-_tmp54 * _tmp59 * _tmp95 + _tmp56 * _tmp57 * _tmp95);
  const Scalar _tmp97 = _tmp72 * (_tmp63 * _tmp70 - _tmp66 * _tmp71 + _tmp70 * _tmp96);
  const Scalar _tmp98 = -_tmp42 * _tmp53 + _tmp48 * _tmp52 + _tmp52 * _tmp96 - _tmp62 * _tmp97;
  const Scalar _tmp99 = Scalar(1.0) / (_tmp98);
  const Scalar _tmp100 = Scalar(1.0) * _tmp99;
  const Scalar _tmp101 = _tmp31 + Scalar(-4.7744369927459998);
  const Scalar _tmp102 = _tmp20 + Scalar(-2.7171519410699099);
  const Scalar _tmp103 =
      std::pow(Scalar(std::pow(_tmp101, Scalar(2)) + std::pow(_tmp102, Scalar(2))),
               Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp104 = _tmp101 * _tmp103;
  const Scalar _tmp105 = _tmp102 * _tmp103;
  const Scalar _tmp106 = fh1 * (_tmp104 * _tmp19 - _tmp105 * _tmp30);
  const Scalar _tmp107 = Scalar(1.0) * _tmp77;
  const Scalar _tmp108 = _tmp107 * _tmp84 - Scalar(1.0) * _tmp85;
  const Scalar _tmp109 = _tmp87 * _tmp98;
  const Scalar _tmp110 = _tmp86 * _tmp99;
  const Scalar _tmp111 = _tmp110 * (-_tmp108 * _tmp109 - Scalar(1.0) * _tmp97);
  const Scalar _tmp112 = _tmp87 * (_tmp108 + _tmp111);
  const Scalar _tmp113 = -_tmp112 * _tmp62 + Scalar(1.0);
  const Scalar _tmp114 = _tmp104 * fh1;
  const Scalar _tmp115 = _tmp61 * _tmp77 + _tmp78;
  const Scalar _tmp116 = -_tmp115 * _tmp84 + _tmp61 * _tmp85 - _tmp75;
  const Scalar _tmp117 = _tmp110 * (-_tmp109 * _tmp116 + _tmp61 * _tmp97 - _tmp96);
  const Scalar _tmp118 = _tmp87 * (_tmp116 + _tmp117);
  const Scalar _tmp119 = -_tmp118 * _tmp62 - _tmp61;
  const Scalar _tmp120 = _tmp105 * fh1;
  const Scalar _tmp121 = -_tmp106 * _tmp93 * (_tmp100 * _tmp52 - _tmp100 * _tmp91) -
                         _tmp114 * _tmp93 * (_tmp112 * _tmp52 + _tmp113 * _tmp90) -
                         _tmp120 * _tmp93 * (_tmp118 * _tmp52 + _tmp119 * _tmp90 + Scalar(1.0)) -
                         _tmp33 * _tmp94;
  const Scalar _tmp122 = Scalar(1.0) / (_tmp121);
  const Scalar _tmp123 = fh1 * (_tmp34 + _tmp73);
  const Scalar _tmp124 = _tmp105 * _tmp123 + Scalar(5.1796800000000003) * _tmp16 + _tmp19 * fv1;
  const Scalar _tmp125 = _tmp48 + _tmp81;
  const Scalar _tmp126 = _tmp125 * _tmp84;
  const Scalar _tmp127 = Scalar(1.0) / (-_tmp126 - _tmp42 + _tmp83);
  const Scalar _tmp128 = Scalar(1.0) * _tmp127;
  const Scalar _tmp129 = _tmp125 * _tmp128;
  const Scalar _tmp130 = -Scalar(1.0) * _tmp128 + Scalar(1.0) * _tmp129 * _tmp82;
  const Scalar _tmp131 = -_tmp104 * _tmp123 - Scalar(5.1796800000000003) * _tmp26 - _tmp30 * fv1;
  const Scalar _tmp132 = _tmp82 * (_tmp126 * _tmp128 + Scalar(1.0));
  const Scalar _tmp133 = _tmp128 * _tmp84;
  const Scalar _tmp134 = -Scalar(1.0) * _tmp132 + Scalar(1.0) * _tmp133;
  const Scalar _tmp135 = _tmp127 * _tmp88;
  const Scalar _tmp136 = _tmp82 * (-_tmp125 * _tmp135 - _tmp79 * _tmp89 + _tmp81);
  const Scalar _tmp137 = -Scalar(1.0) * _tmp135 - Scalar(1.0) * _tmp136 + Scalar(1.0);
  const Scalar _tmp138 = _tmp125 * _tmp127;
  const Scalar _tmp139 = -_tmp107 + _tmp111 * _tmp138 - _tmp112 * _tmp79;
  const Scalar _tmp140 = Scalar(1.0) * _tmp82;
  const Scalar _tmp141 = Scalar(1.0) * fh1;
  const Scalar _tmp142 = -_tmp100 * _tmp79 + _tmp110 * _tmp129;
  const Scalar _tmp143 = _tmp110 * _tmp128;
  const Scalar _tmp144 = _tmp115 + _tmp117 * _tmp138 - _tmp118 * _tmp79;
  const Scalar _tmp145 = _tmp104 * _tmp141 * (_tmp111 * _tmp128 - _tmp139 * _tmp140) +
                         _tmp105 * _tmp141 * (_tmp117 * _tmp128 - _tmp140 * _tmp144) +
                         Scalar(1.0) * _tmp106 * (-_tmp140 * _tmp142 + _tmp143) +
                         _tmp124 * _tmp130 + _tmp131 * _tmp134 + _tmp137 * _tmp33;
  const Scalar _tmp146 = std::asinh(_tmp122 * _tmp145);
  const Scalar _tmp147 = Scalar(1.0) * _tmp146;
  const Scalar _tmp148 = Scalar(9.6622558468725703) * _tmp121;
  const Scalar _tmp149 =
      -_tmp146 * _tmp148 -
      Scalar(4.83288938413423) *
          std::sqrt(
              Scalar(std::pow(Scalar(1 - Scalar(0.20691555724053495) * _tmp55), Scalar(2)) +
                     Scalar(0.13818785160942856) *
                         std::pow(Scalar(-Scalar(0.55661923802822921) * _tmp58 - 1), Scalar(2))));
  const Scalar _tmp150 = Scalar(0.1034955) * _tmp122;
  const Scalar _tmp151 = _tmp149 * _tmp150;
  const Scalar _tmp152 = Scalar(9.6622558468725703) * _tmp94;
  const Scalar _tmp153 = std::pow(_tmp121, Scalar(-2));
  const Scalar _tmp154 = _tmp153 * _tmp94;
  const Scalar _tmp155 = _tmp27 + _tmp45 + _tmp46;
  const Scalar _tmp156 =
      (_tmp122 * (_tmp130 * _tmp19 + _tmp134 * _tmp155 - _tmp137) - _tmp145 * _tmp154) /
      std::sqrt(Scalar(std::pow(_tmp145, Scalar(2)) * _tmp153 + 1));
  const Scalar _tmp157 = _tmp82 * fh1;
  const Scalar _tmp158 = _tmp124 * _tmp128;
  const Scalar _tmp159 = _tmp125 * _tmp82;
  const Scalar _tmp160 = _tmp104 * _tmp139 * _tmp157 + _tmp105 * _tmp144 * _tmp157 +
                         _tmp106 * _tmp142 * _tmp82 + _tmp131 * _tmp132 + _tmp136 * _tmp33 -
                         _tmp158 * _tmp159;
  const Scalar _tmp161 = _tmp72 * fh1;
  const Scalar _tmp162 = _tmp100 * _tmp106;
  const Scalar _tmp163 = _tmp62 * _tmp72;
  const Scalar _tmp164 = _tmp33 * _tmp89;
  const Scalar _tmp165 = _tmp104 * _tmp113 * _tmp161 + _tmp105 * _tmp119 * _tmp161 -
                         _tmp162 * _tmp163 - _tmp163 * _tmp164;
  const Scalar _tmp166 = Scalar(1.0) / (_tmp165);
  const Scalar _tmp167 = std::asinh(_tmp160 * _tmp166);
  const Scalar _tmp168 = Scalar(9.6622558468725703) * _tmp165;
  const Scalar _tmp169 =
      -_tmp167 * _tmp168 -
      Scalar(8.3885017487099702) *
          std::sqrt(
              Scalar(Scalar(0.090199313518583735) *
                         std::pow(Scalar(1 - Scalar(0.39693005512043167) * _tmp67), Scalar(2)) +
                     std::pow(Scalar(-Scalar(0.11921079949155229) * _tmp64 - 1), Scalar(2))));
  const Scalar _tmp170 = Scalar(0.1034955) * _tmp166;
  const Scalar _tmp171 = _tmp169 * _tmp170;
  const Scalar _tmp172 = Scalar(1.0) * _tmp167;
  const Scalar _tmp173 = Scalar(9.6622558468725703) * _tmp89;
  const Scalar _tmp174 = _tmp163 * _tmp173;
  const Scalar _tmp175 = std::pow(_tmp165, Scalar(-2));
  const Scalar _tmp176 = _tmp163 * _tmp175 * _tmp89;
  const Scalar _tmp177 = _tmp128 * _tmp19;
  const Scalar _tmp178 =
      (-_tmp160 * _tmp176 + _tmp166 * (_tmp132 * _tmp155 - _tmp136 - _tmp159 * _tmp177)) /
      std::sqrt(Scalar(std::pow(_tmp160, Scalar(2)) * _tmp175 + 1));
  const Scalar _tmp179 = -_tmp106 * _tmp143 - _tmp111 * _tmp114 * _tmp127 -
                         _tmp117 * _tmp120 * _tmp127 - _tmp131 * _tmp133 + _tmp135 * _tmp33 +
                         _tmp158;
  const Scalar _tmp180 = _tmp112 * _tmp114 + _tmp118 * _tmp120 + _tmp162 + _tmp164;
  const Scalar _tmp181 = Scalar(1.0) / (_tmp180);
  const Scalar _tmp182 = std::asinh(_tmp179 * _tmp181);
  const Scalar _tmp183 = Scalar(1.0) * _tmp182;
  const Scalar _tmp184 = Scalar(9.6622558468725703) * _tmp180;
  const Scalar _tmp185 =
      -_tmp182 * _tmp184 -
      Scalar(8.3196563720703107) *
          std::sqrt(
              Scalar(Scalar(0.057067943527034905) *
                         std::pow(Scalar(-Scalar(0.50315118477277876) * _tmp43 - 1), Scalar(2)) +
                     std::pow(Scalar(-Scalar(0.12019727201198747) * _tmp49 - 1), Scalar(2))));
  const Scalar _tmp186 = Scalar(0.1034955) * _tmp181;
  const Scalar _tmp187 = _tmp185 * _tmp186;
  const Scalar _tmp188 = std::pow(_tmp180, Scalar(-2));
  const Scalar _tmp189 = _tmp188 * _tmp89;
  const Scalar _tmp190 = (_tmp179 * _tmp189 + _tmp181 * (-_tmp133 * _tmp155 - _tmp135 + _tmp177)) /
                         std::sqrt(Scalar(std::pow(_tmp179, Scalar(2)) * _tmp188 + 1));

  // Output terms (1)
  Eigen::Matrix<Scalar, 4, 1> _res;

  _res(0, 0) =
      _tmp32 *
      (-_tmp2 * std::cosh(Scalar(1.0) * _tmp1) +
       _tmp2 *
           std::cosh(
               Scalar(0.1034955) * _tmp0 *
               (-_tmp1 * _tmp32 -
                Scalar(4.7744369927459998) *
                    std::sqrt(Scalar(
                        Scalar(0.32387954179207445) *
                            std::pow(Scalar(1 - Scalar(0.36803241838814449) * _tmp20), Scalar(2)) +
                        std::pow(Scalar(1 - Scalar(0.2094487793051498) * _tmp31), Scalar(2)))))));
  _res(1, 0) = _tmp148 * (-Scalar(1.0) * _tmp156 * std::cosh(_tmp147) -
                          (-Scalar(0.1034955) * _tmp149 * _tmp154 +
                           _tmp150 * (-_tmp146 * _tmp152 - _tmp148 * _tmp156)) *
                              std::cosh(_tmp151)) +
               _tmp152 * (-std::sinh(_tmp147) - std::sinh(_tmp151));
  _res(2, 0) = _tmp168 * (-Scalar(1.0) * _tmp178 * std::cosh(_tmp172) -
                          (-Scalar(0.1034955) * _tmp169 * _tmp176 +
                           _tmp170 * (-_tmp167 * _tmp174 - _tmp168 * _tmp178)) *
                              std::cosh(_tmp171)) +
               _tmp174 * (-std::sinh(_tmp171) - std::sinh(_tmp172));
  _res(3, 0) = -_tmp173 * (-std::sinh(_tmp183) - std::sinh(_tmp187)) +
               _tmp184 * (-Scalar(1.0) * _tmp190 * std::cosh(_tmp183) -
                          (Scalar(0.1034955) * _tmp185 * _tmp189 +
                           _tmp186 * (_tmp173 * _tmp182 - _tmp184 * _tmp190)) *
                              std::cosh(_tmp187));

  return _res;
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym

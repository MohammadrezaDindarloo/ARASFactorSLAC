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
 * Symbolic function: IK_residual_func_cost2_wrt_fv1_Nl0
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
Eigen::Matrix<Scalar, 4, 1> IkResidualFuncCost2WrtFv1Nl0(
    const Scalar fh1, const Scalar fv1, const sym::Rot3<Scalar>& DeltaRot,
    const Eigen::Matrix<Scalar, 3, 1>& position_vector, const sym::Rot3<Scalar>& Rot_init,
    const Scalar epsilon) {
  // Total ops: 596

  // Unused inputs
  (void)epsilon;

  // Input arrays
  const Eigen::Matrix<Scalar, 4, 1>& _DeltaRot = DeltaRot.Data();
  const Eigen::Matrix<Scalar, 4, 1>& _Rot_init = Rot_init.Data();

  // Intermediate terms (192)
  const Scalar _tmp0 = Scalar(1.0) / (fh1);
  const Scalar _tmp1 = std::asinh(_tmp0 * fv1);
  const Scalar _tmp2 = Scalar(1.0) * _tmp0 /
                       std::sqrt(Scalar(1 + std::pow(fv1, Scalar(2)) / std::pow(fh1, Scalar(2))));
  const Scalar _tmp3 = _DeltaRot[0] * _Rot_init[3] - _DeltaRot[1] * _Rot_init[2] +
                       _DeltaRot[2] * _Rot_init[1] + _DeltaRot[3] * _Rot_init[0];
  const Scalar _tmp4 = -2 * std::pow(_tmp3, Scalar(2));
  const Scalar _tmp5 = -_DeltaRot[0] * _Rot_init[1] + _DeltaRot[1] * _Rot_init[0] +
                       _DeltaRot[2] * _Rot_init[3] + _DeltaRot[3] * _Rot_init[2];
  const Scalar _tmp6 = 1 - 2 * std::pow(_tmp5, Scalar(2));
  const Scalar _tmp7 = Scalar(0.20999999999999999) * _tmp4 + Scalar(0.20999999999999999) * _tmp6;
  const Scalar _tmp8 = -_tmp7;
  const Scalar _tmp9 = _DeltaRot[0] * _Rot_init[2] + _DeltaRot[1] * _Rot_init[3] -
                       _DeltaRot[2] * _Rot_init[0] + _DeltaRot[3] * _Rot_init[1];
  const Scalar _tmp10 = 2 * _tmp5;
  const Scalar _tmp11 = _tmp10 * _tmp9;
  const Scalar _tmp12 = -2 * _DeltaRot[0] * _Rot_init[0] - 2 * _DeltaRot[1] * _Rot_init[1] -
                        2 * _DeltaRot[2] * _Rot_init[2] + 2 * _DeltaRot[3] * _Rot_init[3];
  const Scalar _tmp13 = _tmp12 * _tmp3;
  const Scalar _tmp14 = _tmp11 - _tmp13;
  const Scalar _tmp15 = Scalar(0.010999999999999999) * _tmp14;
  const Scalar _tmp16 = -_tmp15;
  const Scalar _tmp17 = 2 * _tmp3 * _tmp9;
  const Scalar _tmp18 = _tmp12 * _tmp5;
  const Scalar _tmp19 = Scalar(0.20999999999999999) * _tmp17 + Scalar(0.20999999999999999) * _tmp18;
  const Scalar _tmp20 = _tmp16 - _tmp19;
  const Scalar _tmp21 = _tmp20 + _tmp8;
  const Scalar _tmp22 = _tmp21 + position_vector(1, 0);
  const Scalar _tmp23 = Scalar(0.20999999999999999) * _tmp17 - Scalar(0.20999999999999999) * _tmp18;
  const Scalar _tmp24 = -_tmp23;
  const Scalar _tmp25 = _tmp10 * _tmp3;
  const Scalar _tmp26 = _tmp12 * _tmp9;
  const Scalar _tmp27 = _tmp25 + _tmp26;
  const Scalar _tmp28 = -Scalar(0.010999999999999999) * _tmp27;
  const Scalar _tmp29 = -2 * std::pow(_tmp9, Scalar(2));
  const Scalar _tmp30 = Scalar(0.20999999999999999) * _tmp29 + Scalar(0.20999999999999999) * _tmp6;
  const Scalar _tmp31 = _tmp28 - _tmp30;
  const Scalar _tmp32 = _tmp24 + _tmp31;
  const Scalar _tmp33 = _tmp32 + position_vector(0, 0);
  const Scalar _tmp34 = Scalar(9.6622558468725703) * fh1;
  const Scalar _tmp35 = Scalar(43.164000000000001) - fv1;
  const Scalar _tmp36 = _tmp23 + _tmp31;
  const Scalar _tmp37 = _tmp36 + position_vector(0, 0);
  const Scalar _tmp38 = _tmp37 + Scalar(1.7965602546229);
  const Scalar _tmp39 = _tmp20 + _tmp7;
  const Scalar _tmp40 = _tmp39 + position_vector(1, 0);
  const Scalar _tmp41 = _tmp40 + Scalar(-4.83288938413423);
  const Scalar _tmp42 = std::pow(Scalar(std::pow(_tmp38, Scalar(2)) + std::pow(_tmp41, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp43 = _tmp38 * _tmp42;
  const Scalar _tmp44 = _tmp19 + _tmp7;
  const Scalar _tmp45 = _tmp16 + _tmp44;
  const Scalar _tmp46 = _tmp45 + position_vector(1, 0);
  const Scalar _tmp47 = _tmp46 + Scalar(-4.7744369927459998);
  const Scalar _tmp48 = _tmp28 + _tmp30;
  const Scalar _tmp49 = _tmp23 + _tmp48;
  const Scalar _tmp50 = _tmp49 + position_vector(0, 0);
  const Scalar _tmp51 = _tmp50 + Scalar(-2.7171519410699099);
  const Scalar _tmp52 = std::pow(Scalar(std::pow(_tmp47, Scalar(2)) + std::pow(_tmp51, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp53 = _tmp51 * _tmp52;
  const Scalar _tmp54 = Scalar(0.20999999999999999) * _tmp11 + Scalar(0.20999999999999999) * _tmp13;
  const Scalar _tmp55 = -_tmp54;
  const Scalar _tmp56 = -Scalar(0.010999999999999999) * _tmp29 -
                        Scalar(0.010999999999999999) * _tmp4 + Scalar(-0.010999999999999999);
  const Scalar _tmp57 = Scalar(0.20999999999999999) * _tmp25 - Scalar(0.20999999999999999) * _tmp26;
  const Scalar _tmp58 = _tmp56 + _tmp57;
  const Scalar _tmp59 = _tmp55 + _tmp58;
  const Scalar _tmp60 = _tmp24 + _tmp48;
  const Scalar _tmp61 = _tmp60 + position_vector(0, 0);
  const Scalar _tmp62 = _tmp61 + Scalar(-2.5193355532036801);
  const Scalar _tmp63 = Scalar(1.0) / (_tmp62);
  const Scalar _tmp64 = _tmp16 + _tmp19 + _tmp8;
  const Scalar _tmp65 = _tmp64 + position_vector(1, 0);
  const Scalar _tmp66 = _tmp65 + Scalar(8.3885017487099702);
  const Scalar _tmp67 = _tmp63 * _tmp66;
  const Scalar _tmp68 = _tmp59 * _tmp67;
  const Scalar _tmp69 = _tmp54 + _tmp58;
  const Scalar _tmp70 = _tmp47 * _tmp52;
  const Scalar _tmp71 = -_tmp53 * _tmp68 + _tmp69 * _tmp70;
  const Scalar _tmp72 = Scalar(1.0) / (_tmp53 * _tmp67 - _tmp70);
  const Scalar _tmp73 = _tmp41 * _tmp42;
  const Scalar _tmp74 = _tmp43 * _tmp67 - _tmp73;
  const Scalar _tmp75 = _tmp72 * _tmp74;
  const Scalar _tmp76 = _tmp43 * _tmp59;
  const Scalar _tmp77 = _tmp56 - _tmp57;
  const Scalar _tmp78 = _tmp54 + _tmp77;
  const Scalar _tmp79 = -_tmp67 * _tmp76 - _tmp71 * _tmp75 + _tmp73 * _tmp78;
  const Scalar _tmp80 = Scalar(1.0) * _tmp64;
  const Scalar _tmp81 = -_tmp80;
  const Scalar _tmp82 = Scalar(1.0) / (_tmp45 + _tmp81);
  const Scalar _tmp83 = Scalar(1.0) * _tmp60;
  const Scalar _tmp84 = _tmp82 * (-_tmp49 + _tmp83);
  const Scalar _tmp85 = _tmp53 * _tmp59 - _tmp53 * _tmp69;
  const Scalar _tmp86 = -_tmp43 * _tmp78 - _tmp75 * _tmp85 + _tmp76 - _tmp79 * _tmp84;
  const Scalar _tmp87 = Scalar(1.0) / (_tmp86);
  const Scalar _tmp88 = _tmp80 * _tmp84 + _tmp83;
  const Scalar _tmp89 = 0;
  const Scalar _tmp90 = _tmp87 * _tmp89;
  const Scalar _tmp91 = _tmp75 * _tmp90;
  const Scalar _tmp92 =
      std::sqrt(Scalar(std::pow(_tmp62, Scalar(2)) + std::pow(_tmp66, Scalar(2))));
  const Scalar _tmp93 = _tmp63 * _tmp92;
  const Scalar _tmp94 = _tmp93 * (_tmp43 * _tmp90 - _tmp53 * _tmp91);
  const Scalar _tmp95 = Scalar(1.0) / (_tmp92);
  const Scalar _tmp96 = _tmp93 * (_tmp60 * _tmp66 * _tmp95 - _tmp62 * _tmp64 * _tmp95);
  const Scalar _tmp97 = _tmp45 * _tmp53 - _tmp49 * _tmp70 + _tmp53 * _tmp96;
  const Scalar _tmp98 = Scalar(1.0) * _tmp72;
  const Scalar _tmp99 = _tmp71 * _tmp98;
  const Scalar _tmp100 = _tmp84 * _tmp99 - _tmp85 * _tmp98;
  const Scalar _tmp101 = -_tmp36 * _tmp73 + _tmp39 * _tmp43 + _tmp43 * _tmp96 - _tmp75 * _tmp97;
  const Scalar _tmp102 = _tmp101 * _tmp87;
  const Scalar _tmp103 = Scalar(1.0) / (_tmp101);
  const Scalar _tmp104 = _tmp103 * _tmp86;
  const Scalar _tmp105 = _tmp104 * (-_tmp100 * _tmp102 - _tmp97 * _tmp98);
  const Scalar _tmp106 = _tmp100 + _tmp105;
  const Scalar _tmp107 = _tmp43 * _tmp87;
  const Scalar _tmp108 = _tmp74 * _tmp87;
  const Scalar _tmp109 = -_tmp106 * _tmp108 + Scalar(1.0);
  const Scalar _tmp110 = _tmp53 * _tmp72;
  const Scalar _tmp111 = _tmp33 + Scalar(1.9874742031097401);
  const Scalar _tmp112 = _tmp22 + Scalar(8.3196563720703107);
  const Scalar _tmp113 =
      std::pow(Scalar(std::pow(_tmp111, Scalar(2)) + std::pow(_tmp112, Scalar(2))),
               Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp114 = _tmp112 * _tmp113;
  const Scalar _tmp115 = _tmp114 * fh1;
  const Scalar _tmp116 = Scalar(1.0) * _tmp103;
  const Scalar _tmp117 = _tmp111 * _tmp113;
  const Scalar _tmp118 = fh1 * (_tmp114 * _tmp32 - _tmp117 * _tmp21);
  const Scalar _tmp119 = _tmp67 * _tmp72;
  const Scalar _tmp120 = _tmp119 * _tmp71 + _tmp68;
  const Scalar _tmp121 = _tmp119 * _tmp85 - _tmp120 * _tmp84 - _tmp59;
  const Scalar _tmp122 = _tmp104 * (-_tmp102 * _tmp121 + _tmp119 * _tmp97 - _tmp96);
  const Scalar _tmp123 = _tmp121 + _tmp122;
  const Scalar _tmp124 = -_tmp108 * _tmp123 - _tmp67;
  const Scalar _tmp125 = _tmp117 * fh1;
  const Scalar _tmp126 = -_tmp115 * _tmp93 * (_tmp106 * _tmp107 + _tmp109 * _tmp110) -
                         _tmp118 * _tmp93 * (_tmp116 * _tmp43 - _tmp116 * _tmp53 * _tmp75) -
                         _tmp125 * _tmp93 * (_tmp107 * _tmp123 + _tmp110 * _tmp124 + Scalar(1.0)) -
                         _tmp35 * _tmp94;
  const Scalar _tmp127 = Scalar(1.0) / (_tmp126);
  const Scalar _tmp128 = fh1 * (_tmp55 + _tmp77);
  const Scalar _tmp129 = _tmp117 * _tmp128 + Scalar(5.1796800000000003) * _tmp27 + _tmp32 * fv1;
  const Scalar _tmp130 = _tmp39 + _tmp81;
  const Scalar _tmp131 = _tmp130 * _tmp84;
  const Scalar _tmp132 = Scalar(1.0) / (-_tmp131 - _tmp36 + _tmp83);
  const Scalar _tmp133 = Scalar(1.0) * _tmp132;
  const Scalar _tmp134 = _tmp130 * _tmp82;
  const Scalar _tmp135 = Scalar(1.0) * _tmp133 * _tmp134 - Scalar(1.0) * _tmp133;
  const Scalar _tmp136 = _tmp130 * _tmp132;
  const Scalar _tmp137 = _tmp79 * _tmp87;
  const Scalar _tmp138 = _tmp120 + _tmp122 * _tmp136 - _tmp123 * _tmp137;
  const Scalar _tmp139 = Scalar(1.0) * _tmp82;
  const Scalar _tmp140 = _tmp132 * _tmp88;
  const Scalar _tmp141 = _tmp82 * (-_tmp130 * _tmp140 - _tmp137 * _tmp89 + _tmp81);
  const Scalar _tmp142 = -Scalar(1.0) * _tmp140 - Scalar(1.0) * _tmp141 + Scalar(1.0);
  const Scalar _tmp143 = _tmp105 * _tmp136 - _tmp106 * _tmp137 - _tmp99;
  const Scalar _tmp144 = _tmp104 * _tmp133;
  const Scalar _tmp145 = -_tmp116 * _tmp79 + _tmp130 * _tmp144;
  const Scalar _tmp146 = -_tmp114 * _tmp128 - Scalar(5.1796800000000003) * _tmp14 - _tmp21 * fv1;
  const Scalar _tmp147 = _tmp82 * (_tmp131 * _tmp133 + Scalar(1.0));
  const Scalar _tmp148 = _tmp133 * _tmp84;
  const Scalar _tmp149 = -Scalar(1.0) * _tmp147 + Scalar(1.0) * _tmp148;
  const Scalar _tmp150 = Scalar(1.0) * _tmp115 * (_tmp105 * _tmp133 - _tmp139 * _tmp143) +
                         Scalar(1.0) * _tmp118 * (-_tmp139 * _tmp145 + _tmp144) +
                         Scalar(1.0) * _tmp125 * (_tmp122 * _tmp133 - _tmp138 * _tmp139) +
                         _tmp129 * _tmp135 + _tmp142 * _tmp35 + _tmp146 * _tmp149;
  const Scalar _tmp151 = std::asinh(_tmp127 * _tmp150);
  const Scalar _tmp152 = Scalar(1.0) * _tmp151;
  const Scalar _tmp153 = Scalar(9.6622558468725703) * _tmp126;
  const Scalar _tmp154 =
      -_tmp151 * _tmp153 -
      Scalar(8.3885017487099702) *
          std::sqrt(
              Scalar(Scalar(0.090199313518583735) *
                         std::pow(Scalar(1 - Scalar(0.39693005512043167) * _tmp61), Scalar(2)) +
                     std::pow(Scalar(-Scalar(0.11921079949155229) * _tmp65 - 1), Scalar(2))));
  const Scalar _tmp155 = Scalar(0.1034955) * _tmp127;
  const Scalar _tmp156 = _tmp154 * _tmp155;
  const Scalar _tmp157 = Scalar(9.6622558468725703) * _tmp94;
  const Scalar _tmp158 = std::pow(_tmp126, Scalar(-2));
  const Scalar _tmp159 = _tmp158 * _tmp94;
  const Scalar _tmp160 = _tmp15 + _tmp44;
  const Scalar _tmp161 =
      (_tmp127 * (_tmp135 * _tmp32 - _tmp142 + _tmp149 * _tmp160) - _tmp150 * _tmp159) /
      std::sqrt(Scalar(std::pow(_tmp150, Scalar(2)) * _tmp158 + 1));
  const Scalar _tmp162 = _tmp35 * _tmp90;
  const Scalar _tmp163 = _tmp116 * _tmp118;
  const Scalar _tmp164 =
      _tmp109 * _tmp115 * _tmp72 + _tmp124 * _tmp125 * _tmp72 - _tmp162 * _tmp75 - _tmp163 * _tmp75;
  const Scalar _tmp165 = Scalar(1.0) / (_tmp164);
  const Scalar _tmp166 = _tmp129 * _tmp133;
  const Scalar _tmp167 = _tmp115 * _tmp143 * _tmp82 + _tmp118 * _tmp145 * _tmp82 +
                         _tmp125 * _tmp138 * _tmp82 - _tmp134 * _tmp166 + _tmp141 * _tmp35 +
                         _tmp146 * _tmp147;
  const Scalar _tmp168 = std::asinh(_tmp165 * _tmp167);
  const Scalar _tmp169 = Scalar(1.0) * _tmp168;
  const Scalar _tmp170 = Scalar(9.6622558468725703) * _tmp164;
  const Scalar _tmp171 =
      -_tmp168 * _tmp170 -
      Scalar(4.7744369927459998) *
          std::sqrt(
              Scalar(std::pow(Scalar(1 - Scalar(0.2094487793051498) * _tmp46), Scalar(2)) +
                     Scalar(0.32387954179207445) *
                         std::pow(Scalar(1 - Scalar(0.36803241838814449) * _tmp50), Scalar(2))));
  const Scalar _tmp172 = Scalar(0.1034955) * _tmp165;
  const Scalar _tmp173 = _tmp171 * _tmp172;
  const Scalar _tmp174 = Scalar(9.6622558468725703) * _tmp90;
  const Scalar _tmp175 = _tmp174 * _tmp75;
  const Scalar _tmp176 = std::pow(_tmp164, Scalar(-2));
  const Scalar _tmp177 = _tmp176 * _tmp91;
  const Scalar _tmp178 = _tmp133 * _tmp32;
  const Scalar _tmp179 =
      (_tmp165 * (-_tmp134 * _tmp178 - _tmp141 + _tmp147 * _tmp160) - _tmp167 * _tmp177) /
      std::sqrt(Scalar(std::pow(_tmp167, Scalar(2)) * _tmp176 + 1));
  const Scalar _tmp180 =
      _tmp106 * _tmp115 * _tmp87 + _tmp123 * _tmp125 * _tmp87 + _tmp162 + _tmp163;
  const Scalar _tmp181 = Scalar(1.0) / (_tmp180);
  const Scalar _tmp182 = -_tmp105 * _tmp115 * _tmp132 - _tmp118 * _tmp144 -
                         _tmp122 * _tmp125 * _tmp132 + _tmp140 * _tmp35 - _tmp146 * _tmp148 +
                         _tmp166;
  const Scalar _tmp183 = std::asinh(_tmp181 * _tmp182);
  const Scalar _tmp184 = Scalar(1.0) * _tmp183;
  const Scalar _tmp185 = std::pow(_tmp180, Scalar(-2));
  const Scalar _tmp186 = _tmp185 * _tmp90;
  const Scalar _tmp187 = (_tmp181 * (-_tmp140 - _tmp148 * _tmp160 + _tmp178) + _tmp182 * _tmp186) /
                         std::sqrt(Scalar(std::pow(_tmp182, Scalar(2)) * _tmp185 + 1));
  const Scalar _tmp188 = Scalar(9.6622558468725703) * _tmp180;
  const Scalar _tmp189 =
      -_tmp183 * _tmp188 -
      Scalar(4.83288938413423) *
          std::sqrt(
              Scalar(std::pow(Scalar(1 - Scalar(0.20691555724053495) * _tmp40), Scalar(2)) +
                     Scalar(0.13818785160942856) *
                         std::pow(Scalar(-Scalar(0.55661923802822921) * _tmp37 - 1), Scalar(2))));
  const Scalar _tmp190 = Scalar(0.1034955) * _tmp181;
  const Scalar _tmp191 = _tmp189 * _tmp190;

  // Output terms (1)
  Eigen::Matrix<Scalar, 4, 1> _res;

  _res(0, 0) =
      _tmp34 *
      (-_tmp2 * std::cosh(Scalar(1.0) * _tmp1) +
       _tmp2 * std::cosh(
                   Scalar(0.1034955) * _tmp0 *
                   (-_tmp1 * _tmp34 -
                    Scalar(8.3196563720703107) *
                        std::sqrt(Scalar(
                            std::pow(Scalar(-Scalar(0.12019727201198747) * _tmp22 - 1), Scalar(2)) +
                            Scalar(0.057067943527034905) *
                                std::pow(Scalar(-Scalar(0.50315118477277876) * _tmp33 - 1),
                                         Scalar(2)))))));
  _res(1, 0) = _tmp153 * (-Scalar(1.0) * _tmp161 * std::cosh(_tmp152) -
                          (-Scalar(0.1034955) * _tmp154 * _tmp159 +
                           _tmp155 * (-_tmp151 * _tmp157 - _tmp153 * _tmp161)) *
                              std::cosh(_tmp156)) +
               _tmp157 * (-std::sinh(_tmp152) - std::sinh(_tmp156));
  _res(2, 0) = _tmp170 * (-Scalar(1.0) * _tmp179 * std::cosh(_tmp169) -
                          (-Scalar(0.1034955) * _tmp171 * _tmp177 +
                           _tmp172 * (-_tmp168 * _tmp175 - _tmp170 * _tmp179)) *
                              std::cosh(_tmp173)) +
               _tmp175 * (-std::sinh(_tmp169) - std::sinh(_tmp173));
  _res(3, 0) = -_tmp174 * (-std::sinh(_tmp184) - std::sinh(_tmp191)) +
               _tmp188 * (-Scalar(1.0) * _tmp187 * std::cosh(_tmp184) -
                          (Scalar(0.1034955) * _tmp186 * _tmp189 +
                           _tmp190 * (_tmp174 * _tmp183 - _tmp187 * _tmp188)) *
                              std::cosh(_tmp191));

  return _res;
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym

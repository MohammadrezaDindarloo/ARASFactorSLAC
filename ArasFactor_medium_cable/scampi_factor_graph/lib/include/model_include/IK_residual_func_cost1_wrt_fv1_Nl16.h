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
 * Symbolic function: IK_residual_func_cost1_wrt_fv1_Nl16
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
Eigen::Matrix<Scalar, 4, 1> IkResidualFuncCost1WrtFv1Nl16(
    const Scalar fh1, const Scalar fv1, const sym::Rot3<Scalar>& DeltaRot,
    const Eigen::Matrix<Scalar, 3, 1>& position_vector, const sym::Rot3<Scalar>& Rot_init,
    const Scalar epsilon) {
  // Total ops: 604

  // Unused inputs
  (void)epsilon;

  // Input arrays
  const Eigen::Matrix<Scalar, 4, 1>& _DeltaRot = DeltaRot.Data();
  const Eigen::Matrix<Scalar, 4, 1>& _Rot_init = Rot_init.Data();

  // Intermediate terms (188)
  const Scalar _tmp0 = Scalar(1.0) / (fh1);
  const Scalar _tmp1 = std::asinh(_tmp0 * fv1);
  const Scalar _tmp2 = Scalar(1.0) * _tmp0 /
                       std::sqrt(Scalar(1 + std::pow(fv1, Scalar(2)) / std::pow(fh1, Scalar(2))));
  const Scalar _tmp3 = _DeltaRot[0] * _Rot_init[2] + _DeltaRot[1] * _Rot_init[3] -
                       _DeltaRot[2] * _Rot_init[0] + _DeltaRot[3] * _Rot_init[1];
  const Scalar _tmp4 = _DeltaRot[0] * _Rot_init[3] - _DeltaRot[1] * _Rot_init[2] +
                       _DeltaRot[2] * _Rot_init[1] + _DeltaRot[3] * _Rot_init[0];
  const Scalar _tmp5 = 2 * _tmp3 * _tmp4;
  const Scalar _tmp6 = -_DeltaRot[0] * _Rot_init[1] + _DeltaRot[1] * _Rot_init[0] +
                       _DeltaRot[2] * _Rot_init[3] + _DeltaRot[3] * _Rot_init[2];
  const Scalar _tmp7 = -2 * _DeltaRot[0] * _Rot_init[0] - 2 * _DeltaRot[1] * _Rot_init[1] -
                       2 * _DeltaRot[2] * _Rot_init[2] + 2 * _DeltaRot[3] * _Rot_init[3];
  const Scalar _tmp8 = _tmp6 * _tmp7;
  const Scalar _tmp9 = Scalar(0.20999999999999999) * _tmp5 - Scalar(0.20999999999999999) * _tmp8;
  const Scalar _tmp10 = 2 * _tmp6;
  const Scalar _tmp11 = _tmp10 * _tmp4;
  const Scalar _tmp12 = _tmp3 * _tmp7;
  const Scalar _tmp13 = _tmp11 + _tmp12;
  const Scalar _tmp14 = -Scalar(0.010999999999999999) * _tmp13;
  const Scalar _tmp15 = -2 * std::pow(_tmp6, Scalar(2));
  const Scalar _tmp16 = 1 - 2 * std::pow(_tmp3, Scalar(2));
  const Scalar _tmp17 = Scalar(0.20999999999999999) * _tmp15 + Scalar(0.20999999999999999) * _tmp16;
  const Scalar _tmp18 = _tmp14 + _tmp17;
  const Scalar _tmp19 = _tmp18 + _tmp9;
  const Scalar _tmp20 = _tmp19 + position_vector(0, 0);
  const Scalar _tmp21 = -2 * std::pow(_tmp4, Scalar(2));
  const Scalar _tmp22 = Scalar(0.20999999999999999) * _tmp15 +
                        Scalar(0.20999999999999999) * _tmp21 + Scalar(0.20999999999999999);
  const Scalar _tmp23 = _tmp10 * _tmp3;
  const Scalar _tmp24 = _tmp4 * _tmp7;
  const Scalar _tmp25 = _tmp23 - _tmp24;
  const Scalar _tmp26 = Scalar(0.010999999999999999) * _tmp25;
  const Scalar _tmp27 = -_tmp26;
  const Scalar _tmp28 = Scalar(0.20999999999999999) * _tmp5 + Scalar(0.20999999999999999) * _tmp8;
  const Scalar _tmp29 = _tmp27 + _tmp28;
  const Scalar _tmp30 = _tmp22 + _tmp29;
  const Scalar _tmp31 = _tmp30 + position_vector(1, 0);
  const Scalar _tmp32 = Scalar(9.6622558468725703) * fh1;
  const Scalar _tmp33 = _tmp14 - _tmp17;
  const Scalar _tmp34 = _tmp33 + _tmp9;
  const Scalar _tmp35 = _tmp34 + position_vector(0, 0);
  const Scalar _tmp36 = -_tmp28;
  const Scalar _tmp37 = _tmp27 + _tmp36;
  const Scalar _tmp38 = _tmp22 + _tmp37;
  const Scalar _tmp39 = _tmp38 + position_vector(1, 0);
  const Scalar _tmp40 = Scalar(43.164000000000001) - fv1;
  const Scalar _tmp41 = Scalar(0.20999999999999999) * _tmp23 + Scalar(0.20999999999999999) * _tmp24;
  const Scalar _tmp42 =
      -Scalar(0.010999999999999999) * _tmp16 - Scalar(0.010999999999999999) * _tmp21;
  const Scalar _tmp43 = Scalar(0.20999999999999999) * _tmp11 - Scalar(0.20999999999999999) * _tmp12;
  const Scalar _tmp44 = _tmp42 - _tmp43;
  const Scalar _tmp45 = _tmp41 + _tmp44;
  const Scalar _tmp46 = -_tmp22;
  const Scalar _tmp47 = _tmp29 + _tmp46;
  const Scalar _tmp48 = _tmp47 + position_vector(1, 0);
  const Scalar _tmp49 = _tmp48 + Scalar(8.3885017487099702);
  const Scalar _tmp50 = -_tmp9;
  const Scalar _tmp51 = _tmp18 + _tmp50;
  const Scalar _tmp52 = _tmp51 + position_vector(0, 0);
  const Scalar _tmp53 = _tmp52 + Scalar(-2.5193355532036801);
  const Scalar _tmp54 = std::pow(Scalar(std::pow(_tmp49, Scalar(2)) + std::pow(_tmp53, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp55 = _tmp53 * _tmp54;
  const Scalar _tmp56 = -_tmp41;
  const Scalar _tmp57 = _tmp42 + _tmp43;
  const Scalar _tmp58 = _tmp56 + _tmp57;
  const Scalar _tmp59 = _tmp39 + Scalar(-4.83288938413423);
  const Scalar _tmp60 = _tmp35 + Scalar(1.7965602546229);
  const Scalar _tmp61 = Scalar(1.0) / (_tmp60);
  const Scalar _tmp62 = _tmp59 * _tmp61;
  const Scalar _tmp63 = _tmp49 * _tmp54;
  const Scalar _tmp64 = _tmp55 * _tmp62 - _tmp63;
  const Scalar _tmp65 = _tmp37 + _tmp46;
  const Scalar _tmp66 = _tmp65 + position_vector(1, 0);
  const Scalar _tmp67 = _tmp66 + Scalar(8.3196563720703107);
  const Scalar _tmp68 = _tmp33 + _tmp50;
  const Scalar _tmp69 = _tmp68 + position_vector(0, 0);
  const Scalar _tmp70 = _tmp69 + Scalar(1.9874742031097401);
  const Scalar _tmp71 = std::pow(Scalar(std::pow(_tmp67, Scalar(2)) + std::pow(_tmp70, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp72 = _tmp67 * _tmp71;
  const Scalar _tmp73 = _tmp70 * _tmp71;
  const Scalar _tmp74 = Scalar(1.0) / (_tmp62 * _tmp73 - _tmp72);
  const Scalar _tmp75 = _tmp44 + _tmp56;
  const Scalar _tmp76 = _tmp45 * _tmp73;
  const Scalar _tmp77 = _tmp74 * (-_tmp73 * _tmp75 + _tmp76);
  const Scalar _tmp78 = _tmp45 * _tmp62;
  const Scalar _tmp79 = _tmp74 * (-_tmp62 * _tmp76 + _tmp72 * _tmp75);
  const Scalar _tmp80 = -_tmp55 * _tmp78 + _tmp58 * _tmp63 - _tmp64 * _tmp79;
  const Scalar _tmp81 = Scalar(1.0) * _tmp38;
  const Scalar _tmp82 = -_tmp81;
  const Scalar _tmp83 = Scalar(1.0) / (_tmp65 + _tmp82);
  const Scalar _tmp84 = Scalar(1.0) * _tmp34;
  const Scalar _tmp85 = _tmp83 * (-_tmp68 + _tmp84);
  const Scalar _tmp86 = _tmp45 * _tmp55 - _tmp55 * _tmp58 - _tmp64 * _tmp77 - _tmp80 * _tmp85;
  const Scalar _tmp87 = Scalar(1.0) / (_tmp86);
  const Scalar _tmp88 = _tmp81 * _tmp85 + _tmp84;
  const Scalar _tmp89 = 0;
  const Scalar _tmp90 = _tmp73 * _tmp74;
  const Scalar _tmp91 = _tmp64 * _tmp90;
  const Scalar _tmp92 =
      std::sqrt(Scalar(std::pow(_tmp59, Scalar(2)) + std::pow(_tmp60, Scalar(2))));
  const Scalar _tmp93 = _tmp61 * _tmp92;
  const Scalar _tmp94 = _tmp93 * (_tmp55 * _tmp89 - _tmp89 * _tmp91);
  const Scalar _tmp95 = Scalar(1.0) / (_tmp92);
  const Scalar _tmp96 = _tmp93 * (_tmp34 * _tmp59 * _tmp95 - _tmp38 * _tmp60 * _tmp95);
  const Scalar _tmp97 = _tmp74 * (_tmp65 * _tmp73 - _tmp68 * _tmp72 + _tmp73 * _tmp96);
  const Scalar _tmp98 = Scalar(1.0) * _tmp79;
  const Scalar _tmp99 = -Scalar(1.0) * _tmp77 + _tmp85 * _tmp98;
  const Scalar _tmp100 = _tmp47 * _tmp55 - _tmp51 * _tmp63 + _tmp55 * _tmp96 - _tmp64 * _tmp97;
  const Scalar _tmp101 = _tmp100 * _tmp87;
  const Scalar _tmp102 = Scalar(1.0) / (_tmp100);
  const Scalar _tmp103 = _tmp102 * _tmp86;
  const Scalar _tmp104 = _tmp103 * (-_tmp101 * _tmp99 - Scalar(1.0) * _tmp97);
  const Scalar _tmp105 = _tmp87 * (_tmp104 + _tmp99);
  const Scalar _tmp106 = -_tmp105 * _tmp64 + Scalar(1.0);
  const Scalar _tmp107 = _tmp31 + Scalar(-4.7744369927459998);
  const Scalar _tmp108 = _tmp20 + Scalar(-2.7171519410699099);
  const Scalar _tmp109 =
      std::pow(Scalar(std::pow(_tmp107, Scalar(2)) + std::pow(_tmp108, Scalar(2))),
               Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp110 = _tmp107 * _tmp109;
  const Scalar _tmp111 = _tmp110 * fh1;
  const Scalar _tmp112 = Scalar(1.0) * _tmp102;
  const Scalar _tmp113 = _tmp108 * _tmp109;
  const Scalar _tmp114 = fh1 * (_tmp110 * _tmp19 - _tmp113 * _tmp30);
  const Scalar _tmp115 = _tmp62 * _tmp79 + _tmp78;
  const Scalar _tmp116 = -_tmp115 * _tmp85 - _tmp45 + _tmp62 * _tmp77;
  const Scalar _tmp117 = _tmp103 * (-_tmp101 * _tmp116 + _tmp62 * _tmp97 - _tmp96);
  const Scalar _tmp118 = _tmp87 * (_tmp116 + _tmp117);
  const Scalar _tmp119 = -_tmp118 * _tmp64 - _tmp62;
  const Scalar _tmp120 = _tmp113 * fh1;
  const Scalar _tmp121 = -_tmp111 * _tmp93 * (_tmp105 * _tmp55 + _tmp106 * _tmp90) -
                         _tmp114 * _tmp93 * (_tmp112 * _tmp55 - _tmp112 * _tmp91) -
                         _tmp120 * _tmp93 * (_tmp118 * _tmp55 + _tmp119 * _tmp90 + Scalar(1.0)) -
                         _tmp40 * _tmp94;
  const Scalar _tmp122 = Scalar(1.0) / (_tmp121);
  const Scalar _tmp123 = _tmp47 + _tmp82;
  const Scalar _tmp124 = _tmp123 * _tmp85;
  const Scalar _tmp125 = Scalar(1.0) / (-_tmp124 - _tmp51 + _tmp84);
  const Scalar _tmp126 = _tmp123 * _tmp125;
  const Scalar _tmp127 = _tmp104 * _tmp126 - _tmp105 * _tmp80 - _tmp98;
  const Scalar _tmp128 = Scalar(1.0) * _tmp83;
  const Scalar _tmp129 = Scalar(1.0) * _tmp125;
  const Scalar _tmp130 = fh1 * (_tmp41 + _tmp57);
  const Scalar _tmp131 = -_tmp110 * _tmp130 - Scalar(5.1796800000000003) * _tmp25 - _tmp30 * fv1;
  const Scalar _tmp132 = _tmp83 * (_tmp124 * _tmp129 + Scalar(1.0));
  const Scalar _tmp133 = _tmp129 * _tmp85;
  const Scalar _tmp134 = -Scalar(1.0) * _tmp132 + Scalar(1.0) * _tmp133;
  const Scalar _tmp135 = _tmp123 * _tmp129;
  const Scalar _tmp136 = _tmp103 * _tmp135 - _tmp112 * _tmp80;
  const Scalar _tmp137 = _tmp103 * _tmp129;
  const Scalar _tmp138 = _tmp125 * _tmp88;
  const Scalar _tmp139 = _tmp83 * (-_tmp123 * _tmp138 - _tmp80 * _tmp89 + _tmp82);
  const Scalar _tmp140 = -Scalar(1.0) * _tmp138 - Scalar(1.0) * _tmp139 + Scalar(1.0);
  const Scalar _tmp141 = _tmp115 + _tmp117 * _tmp126 - _tmp118 * _tmp80;
  const Scalar _tmp142 = _tmp113 * _tmp130 + Scalar(5.1796800000000003) * _tmp13 + _tmp19 * fv1;
  const Scalar _tmp143 = -Scalar(1.0) * _tmp129 + Scalar(1.0) * _tmp135 * _tmp83;
  const Scalar _tmp144 = Scalar(1.0) * _tmp111 * (_tmp104 * _tmp129 - _tmp127 * _tmp128) +
                         Scalar(1.0) * _tmp114 * (-_tmp128 * _tmp136 + _tmp137) +
                         Scalar(1.0) * _tmp120 * (_tmp117 * _tmp129 - _tmp128 * _tmp141) +
                         _tmp131 * _tmp134 + _tmp140 * _tmp40 + _tmp142 * _tmp143;
  const Scalar _tmp145 = std::asinh(_tmp122 * _tmp144);
  const Scalar _tmp146 = Scalar(9.6622558468725703) * _tmp121;
  const Scalar _tmp147 =
      -_tmp145 * _tmp146 -
      Scalar(4.83288938413423) *
          std::sqrt(
              Scalar(std::pow(Scalar(1 - Scalar(0.20691555724053495) * _tmp39), Scalar(2)) +
                     Scalar(0.13818785160942856) *
                         std::pow(Scalar(-Scalar(0.55661923802822921) * _tmp35 - 1), Scalar(2))));
  const Scalar _tmp148 = Scalar(0.1034955) * _tmp122;
  const Scalar _tmp149 = _tmp147 * _tmp148;
  const Scalar _tmp150 = Scalar(1.0) * _tmp145;
  const Scalar _tmp151 = Scalar(9.6622558468725703) * _tmp94;
  const Scalar _tmp152 = std::pow(_tmp121, Scalar(-2));
  const Scalar _tmp153 = _tmp152 * _tmp94;
  const Scalar _tmp154 = _tmp26 + _tmp36 + _tmp46;
  const Scalar _tmp155 =
      (_tmp122 * (_tmp134 * _tmp154 - _tmp140 + _tmp143 * _tmp19) - _tmp144 * _tmp153) /
      std::sqrt(Scalar(std::pow(_tmp144, Scalar(2)) * _tmp152 + 1));
  const Scalar _tmp156 = _tmp112 * _tmp114;
  const Scalar _tmp157 = _tmp64 * _tmp74;
  const Scalar _tmp158 = _tmp40 * _tmp89;
  const Scalar _tmp159 = _tmp106 * _tmp111 * _tmp74 + _tmp119 * _tmp120 * _tmp74 -
                         _tmp156 * _tmp157 - _tmp157 * _tmp158;
  const Scalar _tmp160 = Scalar(1.0) / (_tmp159);
  const Scalar _tmp161 = _tmp129 * _tmp142;
  const Scalar _tmp162 = _tmp123 * _tmp83;
  const Scalar _tmp163 = _tmp111 * _tmp127 * _tmp83 + _tmp114 * _tmp136 * _tmp83 +
                         _tmp120 * _tmp141 * _tmp83 + _tmp131 * _tmp132 + _tmp139 * _tmp40 -
                         _tmp161 * _tmp162;
  const Scalar _tmp164 = std::asinh(_tmp160 * _tmp163);
  const Scalar _tmp165 = Scalar(1.0) * _tmp164;
  const Scalar _tmp166 = Scalar(9.6622558468725703) * _tmp159;
  const Scalar _tmp167 =
      -_tmp164 * _tmp166 -
      Scalar(8.3196563720703107) *
          std::sqrt(
              Scalar(std::pow(Scalar(-Scalar(0.12019727201198747) * _tmp66 - 1), Scalar(2)) +
                     Scalar(0.057067943527034905) *
                         std::pow(Scalar(-Scalar(0.50315118477277876) * _tmp69 - 1), Scalar(2))));
  const Scalar _tmp168 = Scalar(0.1034955) * _tmp160;
  const Scalar _tmp169 = _tmp167 * _tmp168;
  const Scalar _tmp170 = Scalar(9.6622558468725703) * _tmp89;
  const Scalar _tmp171 = _tmp157 * _tmp170;
  const Scalar _tmp172 = _tmp129 * _tmp19;
  const Scalar _tmp173 = std::pow(_tmp159, Scalar(-2));
  const Scalar _tmp174 = _tmp157 * _tmp173 * _tmp89;
  const Scalar _tmp175 =
      (_tmp160 * (_tmp132 * _tmp154 - _tmp139 - _tmp162 * _tmp172) - _tmp163 * _tmp174) /
      std::sqrt(Scalar(std::pow(_tmp163, Scalar(2)) * _tmp173 + 1));
  const Scalar _tmp176 = _tmp105 * _tmp111 + _tmp118 * _tmp120 + _tmp156 + _tmp158;
  const Scalar _tmp177 = Scalar(1.0) / (_tmp176);
  const Scalar _tmp178 = -_tmp104 * _tmp111 * _tmp125 - _tmp114 * _tmp137 -
                         _tmp117 * _tmp120 * _tmp125 - _tmp131 * _tmp133 + _tmp138 * _tmp40 +
                         _tmp161;
  const Scalar _tmp179 = std::asinh(_tmp177 * _tmp178);
  const Scalar _tmp180 = Scalar(9.6622558468725703) * _tmp176;
  const Scalar _tmp181 =
      -_tmp179 * _tmp180 -
      Scalar(8.3885017487099702) *
          std::sqrt(
              Scalar(Scalar(0.090199313518583735) *
                         std::pow(Scalar(1 - Scalar(0.39693005512043167) * _tmp52), Scalar(2)) +
                     std::pow(Scalar(-Scalar(0.11921079949155229) * _tmp48 - 1), Scalar(2))));
  const Scalar _tmp182 = Scalar(0.1034955) * _tmp177;
  const Scalar _tmp183 = _tmp181 * _tmp182;
  const Scalar _tmp184 = Scalar(1.0) * _tmp179;
  const Scalar _tmp185 = std::pow(_tmp176, Scalar(-2));
  const Scalar _tmp186 = _tmp185 * _tmp89;
  const Scalar _tmp187 = (_tmp177 * (-_tmp133 * _tmp154 - _tmp138 + _tmp172) + _tmp178 * _tmp186) /
                         std::sqrt(Scalar(std::pow(_tmp178, Scalar(2)) * _tmp185 + 1));

  // Output terms (1)
  Eigen::Matrix<Scalar, 4, 1> _res;

  _res(0, 0) =
      -_tmp32 *
      (_tmp2 * std::sinh(Scalar(1.0) * _tmp1) +
       _tmp2 *
           std::sinh(
               Scalar(0.1034955) * _tmp0 *
               (-_tmp1 * _tmp32 -
                Scalar(4.7744369927459998) *
                    std::sqrt(Scalar(
                        Scalar(0.32387954179207445) *
                            std::pow(Scalar(1 - Scalar(0.36803241838814449) * _tmp20), Scalar(2)) +
                        std::pow(Scalar(1 - Scalar(0.2094487793051498) * _tmp31), Scalar(2)))))));
  _res(1, 0) =
      -_tmp146 *
          (-Scalar(0.86627065637365697) * _tmp153 + Scalar(1.0) * _tmp155 * std::sinh(_tmp150) -
           (-Scalar(0.1034955) * _tmp147 * _tmp153 +
            _tmp148 * (-_tmp145 * _tmp151 - _tmp146 * _tmp155)) *
               std::sinh(_tmp149)) -
      _tmp151 * (Scalar(0.86627065637365697) * _tmp122 - std::cosh(_tmp149) + std::cosh(_tmp150));
  _res(2, 0) =
      -_tmp166 *
          (-Scalar(0.87679799777269396) * _tmp174 + Scalar(1.0) * _tmp175 * std::sinh(_tmp165) -
           (-Scalar(0.1034955) * _tmp167 * _tmp174 +
            _tmp168 * (-_tmp164 * _tmp171 - _tmp166 * _tmp175)) *
               std::sinh(_tmp169)) -
      _tmp171 * (Scalar(0.87679799777269396) * _tmp160 + std::cosh(_tmp165) - std::cosh(_tmp169));
  _res(3, 0) =
      _tmp170 * (Scalar(0.876505537412406) * _tmp177 - std::cosh(_tmp183) + std::cosh(_tmp184)) -
      _tmp180 * (Scalar(0.876505537412406) * _tmp186 + Scalar(1.0) * _tmp187 * std::sinh(_tmp184) -
                 (Scalar(0.1034955) * _tmp181 * _tmp186 +
                  _tmp182 * (_tmp170 * _tmp179 - _tmp180 * _tmp187)) *
                     std::sinh(_tmp183));

  return _res;
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym

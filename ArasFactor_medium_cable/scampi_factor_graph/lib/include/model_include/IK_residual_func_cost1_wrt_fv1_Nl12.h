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
 * Symbolic function: IK_residual_func_cost1_wrt_fv1_Nl12
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
Eigen::Matrix<Scalar, 4, 1> IkResidualFuncCost1WrtFv1Nl12(
    const Scalar fh1, const Scalar fv1, const sym::Rot3<Scalar>& DeltaRot,
    const Eigen::Matrix<Scalar, 3, 1>& position_vector, const sym::Rot3<Scalar>& Rot_init,
    const Scalar epsilon) {
  // Total ops: 605

  // Unused inputs
  (void)epsilon;

  // Input arrays
  const Eigen::Matrix<Scalar, 4, 1>& _DeltaRot = DeltaRot.Data();
  const Eigen::Matrix<Scalar, 4, 1>& _Rot_init = Rot_init.Data();

  // Intermediate terms (187)
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
  const Scalar _tmp11 = -2 * std::pow(_tmp7, Scalar(2));
  const Scalar _tmp12 = 1 - 2 * std::pow(_tmp4, Scalar(2));
  const Scalar _tmp13 = Scalar(0.20999999999999999) * _tmp11 + Scalar(0.20999999999999999) * _tmp12;
  const Scalar _tmp14 = 2 * _tmp3 * _tmp7;
  const Scalar _tmp15 = _tmp4 * _tmp8;
  const Scalar _tmp16 = _tmp14 + _tmp15;
  const Scalar _tmp17 = -Scalar(0.010999999999999999) * _tmp16;
  const Scalar _tmp18 = _tmp13 + _tmp17;
  const Scalar _tmp19 = _tmp10 + _tmp18;
  const Scalar _tmp20 = _tmp19 + position_vector(0, 0);
  const Scalar _tmp21 = Scalar(0.20999999999999999) * _tmp6 + Scalar(0.20999999999999999) * _tmp9;
  const Scalar _tmp22 = _tmp5 * _tmp7;
  const Scalar _tmp23 = _tmp3 * _tmp8;
  const Scalar _tmp24 = _tmp22 - _tmp23;
  const Scalar _tmp25 = Scalar(0.010999999999999999) * _tmp24;
  const Scalar _tmp26 = -_tmp25;
  const Scalar _tmp27 = -2 * std::pow(_tmp3, Scalar(2));
  const Scalar _tmp28 = Scalar(0.20999999999999999) * _tmp11 +
                        Scalar(0.20999999999999999) * _tmp27 + Scalar(0.20999999999999999);
  const Scalar _tmp29 = _tmp26 + _tmp28;
  const Scalar _tmp30 = _tmp21 + _tmp29;
  const Scalar _tmp31 = _tmp30 + position_vector(1, 0);
  const Scalar _tmp32 = Scalar(9.6622558468725703) * fh1;
  const Scalar _tmp33 = -_tmp13 + _tmp17;
  const Scalar _tmp34 = _tmp10 + _tmp33;
  const Scalar _tmp35 = _tmp34 + position_vector(0, 0);
  const Scalar _tmp36 = _tmp35 + Scalar(1.7965602546229);
  const Scalar _tmp37 = -_tmp21;
  const Scalar _tmp38 = _tmp29 + _tmp37;
  const Scalar _tmp39 = _tmp38 + position_vector(1, 0);
  const Scalar _tmp40 = _tmp39 + Scalar(-4.83288938413423);
  const Scalar _tmp41 = std::pow(Scalar(std::pow(_tmp36, Scalar(2)) + std::pow(_tmp40, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp42 = _tmp36 * _tmp41;
  const Scalar _tmp43 = -_tmp28;
  const Scalar _tmp44 = _tmp37 + _tmp43;
  const Scalar _tmp45 = _tmp26 + _tmp44;
  const Scalar _tmp46 = _tmp45 + position_vector(1, 0);
  const Scalar _tmp47 = _tmp46 + Scalar(8.3196563720703107);
  const Scalar _tmp48 = -_tmp10;
  const Scalar _tmp49 = _tmp33 + _tmp48;
  const Scalar _tmp50 = _tmp49 + position_vector(0, 0);
  const Scalar _tmp51 = _tmp50 + Scalar(1.9874742031097401);
  const Scalar _tmp52 =
      std::sqrt(Scalar(std::pow(_tmp47, Scalar(2)) + std::pow(_tmp51, Scalar(2))));
  const Scalar _tmp53 = Scalar(1.0) / (_tmp52);
  const Scalar _tmp54 = Scalar(1.0) / (_tmp51);
  const Scalar _tmp55 = _tmp52 * _tmp54;
  const Scalar _tmp56 = _tmp55 * (-_tmp45 * _tmp51 * _tmp53 + _tmp47 * _tmp49 * _tmp53);
  const Scalar _tmp57 = _tmp40 * _tmp41;
  const Scalar _tmp58 = _tmp47 * _tmp54;
  const Scalar _tmp59 = _tmp42 * _tmp58 - _tmp57;
  const Scalar _tmp60 = _tmp21 + _tmp26 + _tmp43;
  const Scalar _tmp61 = _tmp60 + position_vector(1, 0);
  const Scalar _tmp62 = _tmp61 + Scalar(8.3885017487099702);
  const Scalar _tmp63 = _tmp18 + _tmp48;
  const Scalar _tmp64 = _tmp63 + position_vector(0, 0);
  const Scalar _tmp65 = _tmp64 + Scalar(-2.5193355532036801);
  const Scalar _tmp66 = std::pow(Scalar(std::pow(_tmp62, Scalar(2)) + std::pow(_tmp65, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp67 = _tmp62 * _tmp66;
  const Scalar _tmp68 = _tmp65 * _tmp66;
  const Scalar _tmp69 = Scalar(1.0) / (_tmp58 * _tmp68 - _tmp67);
  const Scalar _tmp70 = _tmp69 * (_tmp56 * _tmp68 + _tmp60 * _tmp68 - _tmp63 * _tmp67);
  const Scalar _tmp71 = -_tmp34 * _tmp57 + _tmp38 * _tmp42 + _tmp42 * _tmp56 - _tmp59 * _tmp70;
  const Scalar _tmp72 = Scalar(1.0) / (_tmp71);
  const Scalar _tmp73 = Scalar(1.0) * _tmp72;
  const Scalar _tmp74 = _tmp68 * _tmp69;
  const Scalar _tmp75 = _tmp59 * _tmp74;
  const Scalar _tmp76 = _tmp31 + Scalar(-4.7744369927459998);
  const Scalar _tmp77 = _tmp20 + Scalar(-2.7171519410699099);
  const Scalar _tmp78 = std::pow(Scalar(std::pow(_tmp76, Scalar(2)) + std::pow(_tmp77, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp79 = _tmp76 * _tmp78;
  const Scalar _tmp80 = _tmp77 * _tmp78;
  const Scalar _tmp81 = fh1 * (_tmp19 * _tmp79 - _tmp30 * _tmp80);
  const Scalar _tmp82 = Scalar(0.20999999999999999) * _tmp22 + Scalar(0.20999999999999999) * _tmp23;
  const Scalar _tmp83 = -_tmp82;
  const Scalar _tmp84 =
      -Scalar(0.010999999999999999) * _tmp12 - Scalar(0.010999999999999999) * _tmp27;
  const Scalar _tmp85 = Scalar(0.20999999999999999) * _tmp14 - Scalar(0.20999999999999999) * _tmp15;
  const Scalar _tmp86 = _tmp84 - _tmp85;
  const Scalar _tmp87 = _tmp83 + _tmp86;
  const Scalar _tmp88 = _tmp82 + _tmp86;
  const Scalar _tmp89 = _tmp84 + _tmp85;
  const Scalar _tmp90 = _tmp83 + _tmp89;
  const Scalar _tmp91 = _tmp69 * (_tmp68 * _tmp87 - _tmp68 * _tmp90);
  const Scalar _tmp92 = _tmp58 * _tmp87;
  const Scalar _tmp93 = _tmp69 * (_tmp67 * _tmp90 - _tmp68 * _tmp92);
  const Scalar _tmp94 = -_tmp42 * _tmp92 + _tmp57 * _tmp88 - _tmp59 * _tmp93;
  const Scalar _tmp95 = Scalar(1.0) * _tmp45;
  const Scalar _tmp96 = -_tmp95;
  const Scalar _tmp97 = Scalar(1.0) / (_tmp60 + _tmp96);
  const Scalar _tmp98 = Scalar(1.0) * _tmp49;
  const Scalar _tmp99 = -_tmp63 + _tmp98;
  const Scalar _tmp100 = _tmp97 * _tmp99;
  const Scalar _tmp101 = -_tmp100 * _tmp94 + _tmp42 * _tmp87 - _tmp42 * _tmp88 - _tmp59 * _tmp91;
  const Scalar _tmp102 = Scalar(1.0) / (_tmp101);
  const Scalar _tmp103 = Scalar(1.0) * _tmp97;
  const Scalar _tmp104 = _tmp103 * _tmp93 * _tmp99 - Scalar(1.0) * _tmp91;
  const Scalar _tmp105 = _tmp102 * _tmp71;
  const Scalar _tmp106 = _tmp101 * _tmp72;
  const Scalar _tmp107 = _tmp106 * (-_tmp104 * _tmp105 - Scalar(1.0) * _tmp70);
  const Scalar _tmp108 = _tmp102 * (_tmp104 + _tmp107);
  const Scalar _tmp109 = -_tmp108 * _tmp59 + Scalar(1.0);
  const Scalar _tmp110 = _tmp79 * fh1;
  const Scalar _tmp111 = Scalar(43.164000000000001) - fv1;
  const Scalar _tmp112 = _tmp100 * _tmp95 + _tmp98;
  const Scalar _tmp113 = 0;
  const Scalar _tmp114 = _tmp55 * (_tmp113 * _tmp42 - _tmp113 * _tmp75);
  const Scalar _tmp115 = _tmp58 * _tmp93 + _tmp92;
  const Scalar _tmp116 = -_tmp100 * _tmp115 + _tmp58 * _tmp91 - _tmp87;
  const Scalar _tmp117 = _tmp106 * (-_tmp105 * _tmp116 - _tmp56 + _tmp58 * _tmp70);
  const Scalar _tmp118 = _tmp102 * (_tmp116 + _tmp117);
  const Scalar _tmp119 = -_tmp118 * _tmp59 - _tmp58;
  const Scalar _tmp120 = _tmp80 * fh1;
  const Scalar _tmp121 = -_tmp110 * _tmp55 * (_tmp108 * _tmp42 + _tmp109 * _tmp74) -
                         _tmp111 * _tmp114 -
                         _tmp120 * _tmp55 * (_tmp118 * _tmp42 + _tmp119 * _tmp74 + Scalar(1.0)) -
                         _tmp55 * _tmp81 * (_tmp42 * _tmp73 - _tmp73 * _tmp75);
  const Scalar _tmp122 = Scalar(1.0) / (_tmp121);
  const Scalar _tmp123 = fh1 * (_tmp82 + _tmp89);
  const Scalar _tmp124 = _tmp123 * _tmp80 + Scalar(5.1796800000000003) * _tmp16 + _tmp19 * fv1;
  const Scalar _tmp125 = _tmp38 + _tmp96;
  const Scalar _tmp126 = _tmp100 * _tmp125;
  const Scalar _tmp127 = Scalar(1.0) / (-_tmp126 - _tmp34 + _tmp98);
  const Scalar _tmp128 = Scalar(1.0) * _tmp127;
  const Scalar _tmp129 = _tmp125 * _tmp97;
  const Scalar _tmp130 = Scalar(1.0) * _tmp128 * _tmp129 - Scalar(1.0) * _tmp128;
  const Scalar _tmp131 = _tmp106 * _tmp128;
  const Scalar _tmp132 = _tmp125 * _tmp131 - _tmp73 * _tmp94;
  const Scalar _tmp133 = _tmp112 * _tmp127;
  const Scalar _tmp134 = _tmp97 * (-_tmp113 * _tmp94 - _tmp125 * _tmp133 + _tmp96);
  const Scalar _tmp135 = -Scalar(1.0) * _tmp133 - Scalar(1.0) * _tmp134 + Scalar(1.0);
  const Scalar _tmp136 = _tmp125 * _tmp127;
  const Scalar _tmp137 = _tmp115 + _tmp117 * _tmp136 - _tmp118 * _tmp94;
  const Scalar _tmp138 = _tmp107 * _tmp136 - _tmp108 * _tmp94 - Scalar(1.0) * _tmp93;
  const Scalar _tmp139 = -_tmp123 * _tmp79 - Scalar(5.1796800000000003) * _tmp24 - _tmp30 * fv1;
  const Scalar _tmp140 = _tmp100 * _tmp128;
  const Scalar _tmp141 = _tmp126 * _tmp128 + Scalar(1.0);
  const Scalar _tmp142 = -Scalar(1.0) * _tmp103 * _tmp141 + Scalar(1.0) * _tmp140;
  const Scalar _tmp143 =
      Scalar(1.0) * _tmp110 * (-_tmp103 * _tmp138 + _tmp107 * _tmp128) + _tmp111 * _tmp135 +
      Scalar(1.0) * _tmp120 * (-_tmp103 * _tmp137 + _tmp117 * _tmp128) + _tmp124 * _tmp130 +
      _tmp139 * _tmp142 + Scalar(1.0) * _tmp81 * (-_tmp103 * _tmp132 + _tmp131);
  const Scalar _tmp144 = std::asinh(_tmp122 * _tmp143);
  const Scalar _tmp145 = Scalar(9.6622558468725703) * _tmp121;
  const Scalar _tmp146 =
      -_tmp144 * _tmp145 -
      Scalar(8.3196563720703107) *
          std::sqrt(
              Scalar(std::pow(Scalar(-Scalar(0.12019727201198747) * _tmp46 - 1), Scalar(2)) +
                     Scalar(0.057067943527034905) *
                         std::pow(Scalar(-Scalar(0.50315118477277876) * _tmp50 - 1), Scalar(2))));
  const Scalar _tmp147 = Scalar(0.1034955) * _tmp122;
  const Scalar _tmp148 = _tmp146 * _tmp147;
  const Scalar _tmp149 = Scalar(1.0) * _tmp144;
  const Scalar _tmp150 = Scalar(9.6622558468725703) * _tmp114;
  const Scalar _tmp151 = std::pow(_tmp121, Scalar(-2));
  const Scalar _tmp152 = _tmp114 * _tmp151;
  const Scalar _tmp153 = _tmp25 + _tmp44;
  const Scalar _tmp154 =
      (_tmp122 * (_tmp130 * _tmp19 - _tmp135 + _tmp142 * _tmp153) - _tmp143 * _tmp152) /
      std::sqrt(Scalar(std::pow(_tmp143, Scalar(2)) * _tmp151 + 1));
  const Scalar _tmp155 = _tmp111 * _tmp113;
  const Scalar _tmp156 = _tmp59 * _tmp69;
  const Scalar _tmp157 = _tmp73 * _tmp81;
  const Scalar _tmp158 = _tmp109 * _tmp110 * _tmp69 + _tmp119 * _tmp120 * _tmp69 -
                         _tmp155 * _tmp156 - _tmp156 * _tmp157;
  const Scalar _tmp159 = Scalar(1.0) / (_tmp158);
  const Scalar _tmp160 = _tmp141 * _tmp97;
  const Scalar _tmp161 = _tmp124 * _tmp128;
  const Scalar _tmp162 = _tmp110 * _tmp138 * _tmp97 + _tmp111 * _tmp134 +
                         _tmp120 * _tmp137 * _tmp97 - _tmp129 * _tmp161 +
                         _tmp132 * _tmp81 * _tmp97 + _tmp139 * _tmp160;
  const Scalar _tmp163 = std::asinh(_tmp159 * _tmp162);
  const Scalar _tmp164 = Scalar(9.6622558468725703) * _tmp158;
  const Scalar _tmp165 =
      -_tmp163 * _tmp164 -
      Scalar(8.3885017487099702) *
          std::sqrt(
              Scalar(Scalar(0.090199313518583735) *
                         std::pow(Scalar(1 - Scalar(0.39693005512043167) * _tmp64), Scalar(2)) +
                     std::pow(Scalar(-Scalar(0.11921079949155229) * _tmp61 - 1), Scalar(2))));
  const Scalar _tmp166 = Scalar(0.1034955) * _tmp159;
  const Scalar _tmp167 = _tmp165 * _tmp166;
  const Scalar _tmp168 = Scalar(1.0) * _tmp163;
  const Scalar _tmp169 = Scalar(9.6622558468725703) * _tmp113;
  const Scalar _tmp170 = _tmp156 * _tmp169;
  const Scalar _tmp171 = std::pow(_tmp158, Scalar(-2));
  const Scalar _tmp172 = _tmp113 * _tmp156 * _tmp171;
  const Scalar _tmp173 = _tmp128 * _tmp19;
  const Scalar _tmp174 =
      (_tmp159 * (-_tmp129 * _tmp173 - _tmp134 + _tmp153 * _tmp160) - _tmp162 * _tmp172) /
      std::sqrt(Scalar(std::pow(_tmp162, Scalar(2)) * _tmp171 + 1));
  const Scalar _tmp175 = -_tmp107 * _tmp110 * _tmp127 + _tmp111 * _tmp133 -
                         _tmp117 * _tmp120 * _tmp127 - _tmp131 * _tmp81 - _tmp139 * _tmp140 +
                         _tmp161;
  const Scalar _tmp176 = _tmp108 * _tmp110 + _tmp118 * _tmp120 + _tmp155 + _tmp157;
  const Scalar _tmp177 = Scalar(1.0) / (_tmp176);
  const Scalar _tmp178 = std::asinh(_tmp175 * _tmp177);
  const Scalar _tmp179 = Scalar(9.6622558468725703) * _tmp176;
  const Scalar _tmp180 =
      -_tmp178 * _tmp179 -
      Scalar(4.83288938413423) *
          std::sqrt(
              Scalar(std::pow(Scalar(1 - Scalar(0.20691555724053495) * _tmp39), Scalar(2)) +
                     Scalar(0.13818785160942856) *
                         std::pow(Scalar(-Scalar(0.55661923802822921) * _tmp35 - 1), Scalar(2))));
  const Scalar _tmp181 = Scalar(0.1034955) * _tmp177;
  const Scalar _tmp182 = _tmp180 * _tmp181;
  const Scalar _tmp183 = Scalar(1.0) * _tmp178;
  const Scalar _tmp184 = std::pow(_tmp176, Scalar(-2));
  const Scalar _tmp185 = _tmp113 * _tmp184;
  const Scalar _tmp186 = (_tmp175 * _tmp185 + _tmp177 * (-_tmp133 - _tmp140 * _tmp153 + _tmp173)) /
                         std::sqrt(Scalar(std::pow(_tmp175, Scalar(2)) * _tmp184 + 1));

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
      -_tmp145 *
          (-Scalar(0.87679799777269396) * _tmp152 + Scalar(1.0) * _tmp154 * std::sinh(_tmp149) -
           (-Scalar(0.1034955) * _tmp146 * _tmp152 +
            _tmp147 * (-_tmp144 * _tmp150 - _tmp145 * _tmp154)) *
               std::sinh(_tmp148)) -
      _tmp150 * (Scalar(0.87679799777269396) * _tmp122 - std::cosh(_tmp148) + std::cosh(_tmp149));
  _res(2, 0) =
      -_tmp164 *
          (-Scalar(0.876505537412406) * _tmp172 + Scalar(1.0) * _tmp174 * std::sinh(_tmp168) -
           (-Scalar(0.1034955) * _tmp165 * _tmp172 +
            _tmp166 * (-_tmp163 * _tmp170 - _tmp164 * _tmp174)) *
               std::sinh(_tmp167)) -
      _tmp170 * (Scalar(0.876505537412406) * _tmp159 - std::cosh(_tmp167) + std::cosh(_tmp168));
  _res(3, 0) =
      _tmp169 * (Scalar(0.86627065637365697) * _tmp177 - std::cosh(_tmp182) + std::cosh(_tmp183)) -
      _tmp179 *
          (Scalar(0.86627065637365697) * _tmp185 + Scalar(1.0) * _tmp186 * std::sinh(_tmp183) -
           (Scalar(0.1034955) * _tmp180 * _tmp185 +
            _tmp181 * (_tmp169 * _tmp178 - _tmp179 * _tmp186)) *
               std::sinh(_tmp182));

  return _res;
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym

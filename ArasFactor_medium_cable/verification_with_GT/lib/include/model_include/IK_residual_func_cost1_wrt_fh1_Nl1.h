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
 * Symbolic function: IK_residual_func_cost1_wrt_fh1_Nl1
 *
 * Args:
 *     fh1: Scalar
 *     fv1: Scalar
 *     DeltaRot: Rot3
 *     p_init0: Scalar
 *     p_init1: Scalar
 *     p_init2: Scalar
 *     rot_init_x: Scalar
 *     rot_init_y: Scalar
 *     rot_init_z: Scalar
 *     rot_init_w: Scalar
 *     epsilon: Scalar
 *
 * Outputs:
 *     res: Matrix41
 */
template <typename Scalar>
Eigen::Matrix<Scalar, 4, 1> IkResidualFuncCost1WrtFh1Nl1(
    const Scalar fh1, const Scalar fv1, const sym::Rot3<Scalar>& DeltaRot, const Scalar p_init0,
    const Scalar p_init1, const Scalar p_init2, const Scalar rot_init_x, const Scalar rot_init_y,
    const Scalar rot_init_z, const Scalar rot_init_w, const Scalar epsilon) {
  // Total ops: 651

  // Unused inputs
  (void)p_init2;
  (void)epsilon;

  // Input arrays
  const Eigen::Matrix<Scalar, 4, 1>& _DeltaRot = DeltaRot.Data();

  // Intermediate terms (215)
  const Scalar _tmp0 = Scalar(1.0) / (fh1);
  const Scalar _tmp1 = _tmp0 * fv1;
  const Scalar _tmp2 = std::asinh(_tmp1);
  const Scalar _tmp3 = Scalar(1.0) * _tmp2;
  const Scalar _tmp4 = std::pow(fh1, Scalar(-2));
  const Scalar _tmp5 = _DeltaRot[0] * rot_init_z + _DeltaRot[1] * rot_init_w -
                       _DeltaRot[2] * rot_init_x + _DeltaRot[3] * rot_init_y;
  const Scalar _tmp6 = -2 * std::pow(_tmp5, Scalar(2));
  const Scalar _tmp7 = -_DeltaRot[0] * rot_init_y + _DeltaRot[1] * rot_init_x +
                       _DeltaRot[2] * rot_init_w + _DeltaRot[3] * rot_init_z;
  const Scalar _tmp8 = -2 * std::pow(_tmp7, Scalar(2));
  const Scalar _tmp9 = Scalar(0.20999999999999999) * _tmp6 + Scalar(0.20999999999999999) * _tmp8 +
                       Scalar(0.20999999999999999);
  const Scalar _tmp10 = -_tmp9;
  const Scalar _tmp11 = _DeltaRot[0] * rot_init_w - _DeltaRot[1] * rot_init_z +
                        _DeltaRot[2] * rot_init_y + _DeltaRot[3] * rot_init_x;
  const Scalar _tmp12 = 2 * _tmp11;
  const Scalar _tmp13 = _tmp12 * _tmp7;
  const Scalar _tmp14 = -_DeltaRot[0] * rot_init_x - _DeltaRot[1] * rot_init_y -
                        _DeltaRot[2] * rot_init_z + _DeltaRot[3] * rot_init_w;
  const Scalar _tmp15 = 2 * _tmp5;
  const Scalar _tmp16 = _tmp14 * _tmp15;
  const Scalar _tmp17 = _tmp13 + _tmp16;
  const Scalar _tmp18 = -Scalar(0.010999999999999999) * _tmp17;
  const Scalar _tmp19 = _tmp12 * _tmp5;
  const Scalar _tmp20 = 2 * _tmp14 * _tmp7;
  const Scalar _tmp21 = Scalar(0.20999999999999999) * _tmp19 - Scalar(0.20999999999999999) * _tmp20;
  const Scalar _tmp22 = _tmp18 - _tmp21;
  const Scalar _tmp23 = _tmp10 + _tmp22;
  const Scalar _tmp24 = _tmp23 + p_init0;
  const Scalar _tmp25 = Scalar(0.20999999999999999) * _tmp19 + Scalar(0.20999999999999999) * _tmp20;
  const Scalar _tmp26 = -_tmp25;
  const Scalar _tmp27 = _tmp15 * _tmp7;
  const Scalar _tmp28 = _tmp12 * _tmp14;
  const Scalar _tmp29 = _tmp27 - _tmp28;
  const Scalar _tmp30 = -Scalar(0.010999999999999999) * _tmp29;
  const Scalar _tmp31 = 1 - 2 * std::pow(_tmp11, Scalar(2));
  const Scalar _tmp32 = Scalar(0.20999999999999999) * _tmp31 + Scalar(0.20999999999999999) * _tmp8;
  const Scalar _tmp33 = _tmp30 - _tmp32;
  const Scalar _tmp34 = _tmp26 + _tmp33;
  const Scalar _tmp35 = _tmp34 + p_init1;
  const Scalar _tmp36 = Scalar(9.6622558468725703) * _tmp2;
  const Scalar _tmp37 =
      -Scalar(0.1034955) * _tmp36 * fh1 -
      Scalar(0.86104699584133515) *
          std::sqrt(
              Scalar(Scalar(0.057067943376852184) *
                         std::pow(Scalar(-Scalar(0.50315118556004401) * _tmp24 - 1), Scalar(2)) +
                     std::pow(Scalar(-Scalar(0.12019727204189803) * _tmp35 - 1), Scalar(2))));
  const Scalar _tmp38 = _tmp0 * _tmp37;
  const Scalar _tmp39 =
      std::pow(Scalar(_tmp4 * std::pow(fv1, Scalar(2)) + 1), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp40 = _tmp22 + _tmp9;
  const Scalar _tmp41 = Scalar(1.0) * _tmp40;
  const Scalar _tmp42 = _tmp18 + _tmp21;
  const Scalar _tmp43 = _tmp42 + _tmp9;
  const Scalar _tmp44 = _tmp25 + _tmp33;
  const Scalar _tmp45 = Scalar(1.0) * _tmp44;
  const Scalar _tmp46 = -_tmp45;
  const Scalar _tmp47 = _tmp30 + _tmp32;
  const Scalar _tmp48 = _tmp25 + _tmp47;
  const Scalar _tmp49 = _tmp46 + _tmp48;
  const Scalar _tmp50 = _tmp26 + _tmp47;
  const Scalar _tmp51 = Scalar(1.0) / (_tmp46 + _tmp50);
  const Scalar _tmp52 = _tmp10 + _tmp42;
  const Scalar _tmp53 = _tmp51 * (_tmp41 - _tmp52);
  const Scalar _tmp54 = _tmp49 * _tmp53;
  const Scalar _tmp55 = Scalar(1.0) / (_tmp41 - _tmp43 - _tmp54);
  const Scalar _tmp56 = Scalar(1.0) * _tmp55;
  const Scalar _tmp57 = _tmp50 + p_init1;
  const Scalar _tmp58 = _tmp57 + Scalar(-4.8333311099999996);
  const Scalar _tmp59 = _tmp52 + p_init0;
  const Scalar _tmp60 = _tmp59 + Scalar(1.79662371);
  const Scalar _tmp61 = std::pow(Scalar(std::pow(_tmp58, Scalar(2)) + std::pow(_tmp60, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp62 = _tmp58 * _tmp61;
  const Scalar _tmp63 = _tmp60 * _tmp61;
  const Scalar _tmp64 = _tmp40 + p_init0;
  const Scalar _tmp65 = _tmp64 + Scalar(-2.5202214700000001);
  const Scalar _tmp66 = _tmp44 + p_init1;
  const Scalar _tmp67 = _tmp66 + Scalar(8.3888750099999996);
  const Scalar _tmp68 =
      std::sqrt(Scalar(std::pow(_tmp65, Scalar(2)) + std::pow(_tmp67, Scalar(2))));
  const Scalar _tmp69 = Scalar(1.0) / (_tmp68);
  const Scalar _tmp70 = Scalar(1.0) / (_tmp65);
  const Scalar _tmp71 = _tmp68 * _tmp70;
  const Scalar _tmp72 = _tmp71 * (_tmp40 * _tmp67 * _tmp69 - _tmp44 * _tmp65 * _tmp69);
  const Scalar _tmp73 = _tmp50 * _tmp63 - _tmp52 * _tmp62 + _tmp63 * _tmp72;
  const Scalar _tmp74 = _tmp67 * _tmp70;
  const Scalar _tmp75 = Scalar(1.0) / (-_tmp62 + _tmp63 * _tmp74);
  const Scalar _tmp76 = _tmp74 * _tmp75;
  const Scalar _tmp77 = Scalar(0.20999999999999999) * _tmp13 - Scalar(0.20999999999999999) * _tmp16;
  const Scalar _tmp78 =
      -Scalar(0.010999999999999999) * _tmp31 - Scalar(0.010999999999999999) * _tmp6;
  const Scalar _tmp79 = Scalar(0.20999999999999999) * _tmp27 + Scalar(0.20999999999999999) * _tmp28;
  const Scalar _tmp80 = _tmp78 - _tmp79;
  const Scalar _tmp81 = _tmp77 + _tmp80;
  const Scalar _tmp82 = _tmp74 * _tmp81;
  const Scalar _tmp83 = -_tmp77;
  const Scalar _tmp84 = _tmp78 + _tmp79;
  const Scalar _tmp85 = _tmp83 + _tmp84;
  const Scalar _tmp86 = _tmp62 * _tmp85 - _tmp63 * _tmp82;
  const Scalar _tmp87 = _tmp76 * _tmp86 + _tmp82;
  const Scalar _tmp88 = _tmp63 * _tmp81 - _tmp63 * _tmp85;
  const Scalar _tmp89 = -_tmp53 * _tmp87 + _tmp76 * _tmp88 - _tmp81;
  const Scalar _tmp90 = _tmp77 + _tmp84;
  const Scalar _tmp91 = _tmp48 + p_init1;
  const Scalar _tmp92 = _tmp91 + Scalar(-4.7752063900000001);
  const Scalar _tmp93 = _tmp43 + p_init0;
  const Scalar _tmp94 = _tmp93 + Scalar(-2.71799795);
  const Scalar _tmp95 = std::pow(Scalar(std::pow(_tmp92, Scalar(2)) + std::pow(_tmp94, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp96 = _tmp92 * _tmp95;
  const Scalar _tmp97 = _tmp94 * _tmp95;
  const Scalar _tmp98 = _tmp81 * _tmp97;
  const Scalar _tmp99 = _tmp74 * _tmp97 - _tmp96;
  const Scalar _tmp100 = _tmp75 * _tmp99;
  const Scalar _tmp101 = -_tmp100 * _tmp86 - _tmp74 * _tmp98 + _tmp90 * _tmp96;
  const Scalar _tmp102 = -_tmp100 * _tmp88 - _tmp101 * _tmp53 - _tmp90 * _tmp97 + _tmp98;
  const Scalar _tmp103 = Scalar(1.0) / (_tmp102);
  const Scalar _tmp104 = -_tmp100 * _tmp73 - _tmp43 * _tmp96 + _tmp48 * _tmp97 + _tmp72 * _tmp97;
  const Scalar _tmp105 = _tmp103 * _tmp104;
  const Scalar _tmp106 = Scalar(1.0) / (_tmp104);
  const Scalar _tmp107 = _tmp102 * _tmp106;
  const Scalar _tmp108 = _tmp107 * (-_tmp105 * _tmp89 - _tmp72 + _tmp73 * _tmp76);
  const Scalar _tmp109 = _tmp108 + _tmp89;
  const Scalar _tmp110 = _tmp101 * _tmp103;
  const Scalar _tmp111 = _tmp49 * _tmp55;
  const Scalar _tmp112 = _tmp108 * _tmp111 - _tmp109 * _tmp110 + _tmp87;
  const Scalar _tmp113 = Scalar(1.0) * _tmp51;
  const Scalar _tmp114 = _tmp35 + Scalar(8.3196563700000006);
  const Scalar _tmp115 = _tmp24 + Scalar(1.9874742000000001);
  const Scalar _tmp116 =
      std::pow(Scalar(std::pow(_tmp114, Scalar(2)) + std::pow(_tmp115, Scalar(2))),
               Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp117 = _tmp115 * _tmp116;
  const Scalar _tmp118 = Scalar(1.0) * _tmp117 * (_tmp108 * _tmp56 - _tmp112 * _tmp113);
  const Scalar _tmp119 = _tmp80 + _tmp83;
  const Scalar _tmp120 = _tmp119 * fh1;
  const Scalar _tmp121 = _tmp117 * _tmp120 + Scalar(3.29616) * _tmp17 + _tmp23 * fv1;
  const Scalar _tmp122 = _tmp49 * _tmp56;
  const Scalar _tmp123 = Scalar(1.0) * _tmp122 * _tmp51 - Scalar(1.0) * _tmp56;
  const Scalar _tmp124 = Scalar(43.164000000000001) - fv1;
  const Scalar _tmp125 = _tmp41 + _tmp45 * _tmp53;
  const Scalar _tmp126 = _tmp125 * _tmp55;
  const Scalar _tmp127 = 0;
  const Scalar _tmp128 = _tmp51 * (-_tmp110 * _tmp127 - _tmp126 * _tmp49 + _tmp46);
  const Scalar _tmp129 = _tmp114 * _tmp116;
  const Scalar _tmp130 = -_tmp120 * _tmp129 - Scalar(3.29616) * _tmp29 - _tmp34 * fv1;
  const Scalar _tmp131 = _tmp51 * (_tmp54 * _tmp56 + Scalar(1.0));
  const Scalar _tmp132 = _tmp53 * _tmp56;
  const Scalar _tmp133 = -Scalar(1.0) * _tmp131 + Scalar(1.0) * _tmp132;
  const Scalar _tmp134 = Scalar(1.0) * _tmp75;
  const Scalar _tmp135 = _tmp134 * _tmp86;
  const Scalar _tmp136 = -_tmp134 * _tmp88 + _tmp135 * _tmp53;
  const Scalar _tmp137 = _tmp107 * (-_tmp105 * _tmp136 - _tmp134 * _tmp73);
  const Scalar _tmp138 = _tmp136 + _tmp137;
  const Scalar _tmp139 = -_tmp110 * _tmp138 + _tmp111 * _tmp137 - _tmp135;
  const Scalar _tmp140 = Scalar(1.0) * _tmp129 * (-_tmp113 * _tmp139 + _tmp137 * _tmp56);
  const Scalar _tmp141 = -_tmp117 * _tmp34 + _tmp129 * _tmp23;
  const Scalar _tmp142 = _tmp107 * _tmp56;
  const Scalar _tmp143 = Scalar(1.0) * _tmp106;
  const Scalar _tmp144 = -_tmp101 * _tmp143 + _tmp107 * _tmp122;
  const Scalar _tmp145 = Scalar(1.0) * _tmp141 * (-_tmp113 * _tmp144 + _tmp142);
  const Scalar _tmp146 =
      _tmp118 * fh1 + _tmp121 * _tmp123 +
      Scalar(1.0) * _tmp124 * (-_tmp125 * _tmp56 - Scalar(1.0) * _tmp128 + Scalar(1.0)) +
      _tmp130 * _tmp133 + _tmp140 * fh1 + _tmp145 * fh1;
  const Scalar _tmp147 = _tmp103 * _tmp138;
  const Scalar _tmp148 = _tmp103 * _tmp99;
  const Scalar _tmp149 = _tmp75 * (-_tmp138 * _tmp148 + Scalar(1.0));
  const Scalar _tmp150 = _tmp129 * _tmp71 * (_tmp147 * _tmp97 + _tmp149 * _tmp63);
  const Scalar _tmp151 = _tmp103 * _tmp127;
  const Scalar _tmp152 = _tmp100 * _tmp63;
  const Scalar _tmp153 = _tmp103 * _tmp109;
  const Scalar _tmp154 = _tmp75 * (-_tmp109 * _tmp148 - _tmp74);
  const Scalar _tmp155 = _tmp117 * _tmp71 * (_tmp153 * _tmp97 + _tmp154 * _tmp63 + Scalar(1.0));
  const Scalar _tmp156 = _tmp141 * _tmp71 * (-_tmp143 * _tmp152 + _tmp143 * _tmp97);
  const Scalar _tmp157 = -_tmp124 * _tmp71 * (-_tmp151 * _tmp152 + _tmp151 * _tmp97) -
                         _tmp150 * fh1 - _tmp155 * fh1 - _tmp156 * fh1;
  const Scalar _tmp158 = Scalar(1.0) / (_tmp157);
  const Scalar _tmp159 = std::asinh(_tmp146 * _tmp158);
  const Scalar _tmp160 = Scalar(1.0) * _tmp159;
  const Scalar _tmp161 = Scalar(9.6622558468725703) * _tmp157;
  const Scalar _tmp162 =
      -_tmp159 * _tmp161 -
      Scalar(8.3888750099999996) *
          std::sqrt(
              Scalar(Scalar(0.090254729040973036) *
                         std::pow(Scalar(1 - Scalar(0.39679052492160538) * _tmp64), Scalar(2)) +
                     std::pow(Scalar(-Scalar(0.11920549523123722) * _tmp66 - 1), Scalar(2))));
  const Scalar _tmp163 = Scalar(0.1034955) * _tmp158;
  const Scalar _tmp164 = _tmp162 * _tmp163;
  const Scalar _tmp165 = -_tmp150 - _tmp155 - _tmp156;
  const Scalar _tmp166 = Scalar(9.6622558468725703) * _tmp165;
  const Scalar _tmp167 = std::pow(_tmp157, Scalar(-2));
  const Scalar _tmp168 = _tmp165 * _tmp167;
  const Scalar _tmp169 = _tmp117 * _tmp119;
  const Scalar _tmp170 = _tmp119 * _tmp129;
  const Scalar _tmp171 = (-_tmp146 * _tmp168 + _tmp158 * (_tmp118 + _tmp123 * _tmp169 -
                                                          _tmp133 * _tmp170 + _tmp140 + _tmp145)) /
                         std::sqrt(Scalar(std::pow(_tmp146, Scalar(2)) * _tmp167 + 1));
  const Scalar _tmp172 = _tmp141 * _tmp143;
  const Scalar _tmp173 = _tmp172 * fh1;
  const Scalar _tmp174 = _tmp117 * _tmp154;
  const Scalar _tmp175 = _tmp124 * _tmp151;
  const Scalar _tmp176 = _tmp129 * _tmp149;
  const Scalar _tmp177 = -_tmp100 * _tmp173 - _tmp100 * _tmp175 + _tmp174 * fh1 + _tmp176 * fh1;
  const Scalar _tmp178 = std::pow(_tmp177, Scalar(-2));
  const Scalar _tmp179 = -_tmp100 * _tmp172 + _tmp174 + _tmp176;
  const Scalar _tmp180 = _tmp178 * _tmp179;
  const Scalar _tmp181 = Scalar(9.6622558468725703) * _tmp177;
  const Scalar _tmp182 = _tmp129 * _tmp139 * _tmp51;
  const Scalar _tmp183 = _tmp141 * _tmp144 * _tmp51;
  const Scalar _tmp184 = _tmp117 * _tmp51;
  const Scalar _tmp185 = _tmp112 * _tmp184;
  const Scalar _tmp186 = _tmp121 * _tmp56;
  const Scalar _tmp187 = _tmp124 * _tmp128 + _tmp130 * _tmp131 + _tmp182 * fh1 + _tmp183 * fh1 +
                         _tmp185 * fh1 - _tmp186 * _tmp49 * _tmp51;
  const Scalar _tmp188 = Scalar(1.0) / (_tmp177);
  const Scalar _tmp189 =
      (-_tmp180 * _tmp187 +
       _tmp188 * (-_tmp119 * _tmp122 * _tmp184 - _tmp131 * _tmp170 + _tmp182 + _tmp183 + _tmp185)) /
      std::sqrt(Scalar(_tmp178 * std::pow(_tmp187, Scalar(2)) + 1));
  const Scalar _tmp190 = std::asinh(_tmp187 * _tmp188);
  const Scalar _tmp191 = Scalar(9.6622558468725703) * _tmp179;
  const Scalar _tmp192 = Scalar(0.1034955) * _tmp188;
  const Scalar _tmp193 =
      -_tmp181 * _tmp190 -
      Scalar(4.8333311099999996) *
          std::sqrt(
              Scalar(std::pow(Scalar(1 - Scalar(0.20689664689659551) * _tmp57), Scalar(2)) +
                     Scalar(0.13817235445745474) *
                         std::pow(Scalar(-Scalar(0.55659957866191134) * _tmp59 - 1), Scalar(2))));
  const Scalar _tmp194 = _tmp192 * _tmp193;
  const Scalar _tmp195 = Scalar(1.0) * _tmp190;
  const Scalar _tmp196 = _tmp129 * _tmp147;
  const Scalar _tmp197 = _tmp117 * _tmp153;
  const Scalar _tmp198 = _tmp173 + _tmp175 + _tmp196 * fh1 + _tmp197 * fh1;
  const Scalar _tmp199 = std::pow(_tmp198, Scalar(-2));
  const Scalar _tmp200 = _tmp172 + _tmp196 + _tmp197;
  const Scalar _tmp201 = _tmp199 * _tmp200;
  const Scalar _tmp202 = Scalar(1.0) / (_tmp198);
  const Scalar _tmp203 = _tmp141 * _tmp142;
  const Scalar _tmp204 = _tmp108 * _tmp117 * _tmp55;
  const Scalar _tmp205 = _tmp129 * _tmp137 * _tmp55;
  const Scalar _tmp206 = _tmp124 * _tmp126 - _tmp130 * _tmp132 + _tmp186 - _tmp203 * fh1 -
                         _tmp204 * fh1 - _tmp205 * fh1;
  const Scalar _tmp207 = std::asinh(_tmp202 * _tmp206);
  const Scalar _tmp208 = Scalar(1.0) * _tmp207;
  const Scalar _tmp209 = (-_tmp201 * _tmp206 + _tmp202 * (_tmp132 * _tmp170 + _tmp169 * _tmp56 -
                                                          _tmp203 - _tmp204 - _tmp205)) /
                         std::sqrt(Scalar(_tmp199 * std::pow(_tmp206, Scalar(2)) + 1));
  const Scalar _tmp210 = Scalar(9.6622558468725703) * _tmp200;
  const Scalar _tmp211 = Scalar(9.6622558468725703) * _tmp198;
  const Scalar _tmp212 = Scalar(0.1034955) * _tmp202;
  const Scalar _tmp213 =
      -_tmp207 * _tmp211 -
      Scalar(4.7752063900000001) *
          std::sqrt(
              Scalar(std::pow(Scalar(1 - Scalar(0.20941503221602112) * _tmp91), Scalar(2)) +
                     Scalar(0.32397683292140877) *
                         std::pow(Scalar(1 - Scalar(0.36791786395571047) * _tmp93), Scalar(2))));
  const Scalar _tmp214 = _tmp212 * _tmp213;

  // Output terms (1)
  Eigen::Matrix<Scalar, 4, 1> _res;

  _res(0, 0) =
      -Scalar(8.4718465800000011) * _tmp0 -
      Scalar(9.6622558468725703) * fh1 *
          (-Scalar(1.0) * _tmp39 * _tmp4 * fv1 * std::sinh(_tmp3) -
           Scalar(0.87679799772039002) * _tmp4 -
           (Scalar(0.1034955) * _tmp0 * (Scalar(9.6622558468725703) * _tmp1 * _tmp39 - _tmp36) -
            _tmp37 * _tmp4) *
               std::sinh(_tmp38)) -
      Scalar(9.6622558468725703) * std::cosh(_tmp3) +
      Scalar(9.6622558468725703) * std::cosh(_tmp38);
  _res(1, 0) =
      -_tmp161 *
          (-Scalar(0.87653584775870996) * _tmp168 + Scalar(1.0) * _tmp171 * std::sinh(_tmp160) -
           (-Scalar(0.1034955) * _tmp162 * _tmp168 +
            _tmp163 * (-_tmp159 * _tmp166 - _tmp161 * _tmp171)) *
               std::sinh(_tmp164)) -
      _tmp166 * (Scalar(0.87653584775870996) * _tmp158 + std::cosh(_tmp160) - std::cosh(_tmp164));
  _res(2, 0) =
      -_tmp181 *
          (-Scalar(0.86625939559540499) * _tmp180 + Scalar(1.0) * _tmp189 * std::sinh(_tmp195) -
           (-Scalar(0.1034955) * _tmp180 * _tmp193 +
            _tmp192 * (-_tmp181 * _tmp189 - _tmp190 * _tmp191)) *
               std::sinh(_tmp194)) -
      _tmp191 * (Scalar(0.86625939559540499) * _tmp188 - std::cosh(_tmp194) + std::cosh(_tmp195));
  _res(3, 0) =
      -_tmp210 * (Scalar(0.86565325453551001) * _tmp202 + std::cosh(_tmp208) - std::cosh(_tmp214)) -
      _tmp211 *
          (-Scalar(0.86565325453551001) * _tmp201 + Scalar(1.0) * _tmp209 * std::sinh(_tmp208) -
           (-Scalar(0.1034955) * _tmp201 * _tmp213 +
            _tmp212 * (-_tmp207 * _tmp210 - _tmp209 * _tmp211)) *
               std::sinh(_tmp214));

  return _res;
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym

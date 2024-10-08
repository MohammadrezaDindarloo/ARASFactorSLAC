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
 * Symbolic function: IK_residual_func_cost2_wrt_fh1_Nl19
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
Eigen::Matrix<Scalar, 4, 1> IkResidualFuncCost2WrtFh1Nl19(
    const Scalar fh1, const Scalar fv1, const sym::Rot3<Scalar>& DeltaRot, const Scalar p_init0,
    const Scalar p_init1, const Scalar p_init2, const Scalar rot_init_x, const Scalar rot_init_y,
    const Scalar rot_init_z, const Scalar rot_init_w, const Scalar epsilon) {
  // Total ops: 639

  // Unused inputs
  (void)p_init2;
  (void)epsilon;

  // Input arrays
  const Eigen::Matrix<Scalar, 4, 1>& _DeltaRot = DeltaRot.Data();

  // Intermediate terms (210)
  const Scalar _tmp0 = Scalar(1.0) / (fh1);
  const Scalar _tmp1 = _DeltaRot[0] * rot_init_z + _DeltaRot[1] * rot_init_w -
                       _DeltaRot[2] * rot_init_x + _DeltaRot[3] * rot_init_y;
  const Scalar _tmp2 = -2 * std::pow(_tmp1, Scalar(2));
  const Scalar _tmp3 = -_DeltaRot[0] * rot_init_y + _DeltaRot[1] * rot_init_x +
                       _DeltaRot[2] * rot_init_w + _DeltaRot[3] * rot_init_z;
  const Scalar _tmp4 = 1 - 2 * std::pow(_tmp3, Scalar(2));
  const Scalar _tmp5 = Scalar(0.20999999999999999) * _tmp2 + Scalar(0.20999999999999999) * _tmp4;
  const Scalar _tmp6 = -_tmp5;
  const Scalar _tmp7 = _DeltaRot[0] * rot_init_w - _DeltaRot[1] * rot_init_z +
                       _DeltaRot[2] * rot_init_y + _DeltaRot[3] * rot_init_x;
  const Scalar _tmp8 = 2 * _tmp3 * _tmp7;
  const Scalar _tmp9 = -2 * _DeltaRot[0] * rot_init_x - 2 * _DeltaRot[1] * rot_init_y -
                       2 * _DeltaRot[2] * rot_init_z + 2 * _DeltaRot[3] * rot_init_w;
  const Scalar _tmp10 = _tmp1 * _tmp9;
  const Scalar _tmp11 = _tmp10 + _tmp8;
  const Scalar _tmp12 = -Scalar(0.010999999999999999) * _tmp11;
  const Scalar _tmp13 = 2 * _tmp1;
  const Scalar _tmp14 = _tmp13 * _tmp7;
  const Scalar _tmp15 = _tmp3 * _tmp9;
  const Scalar _tmp16 = Scalar(0.20999999999999999) * _tmp14 - Scalar(0.20999999999999999) * _tmp15;
  const Scalar _tmp17 = _tmp12 + _tmp16;
  const Scalar _tmp18 = _tmp17 + _tmp6;
  const Scalar _tmp19 = _tmp18 + p_init0;
  const Scalar _tmp20 = -2 * std::pow(_tmp7, Scalar(2));
  const Scalar _tmp21 = Scalar(0.20999999999999999) * _tmp20 + Scalar(0.20999999999999999) * _tmp4;
  const Scalar _tmp22 = Scalar(0.20999999999999999) * _tmp14 + Scalar(0.20999999999999999) * _tmp15;
  const Scalar _tmp23 = _tmp13 * _tmp3;
  const Scalar _tmp24 = _tmp7 * _tmp9;
  const Scalar _tmp25 = _tmp23 - _tmp24;
  const Scalar _tmp26 = -Scalar(0.010999999999999999) * _tmp25;
  const Scalar _tmp27 = -_tmp22 + _tmp26;
  const Scalar _tmp28 = _tmp21 + _tmp27;
  const Scalar _tmp29 = _tmp28 + p_init1;
  const Scalar _tmp30 = _tmp0 * fv1;
  const Scalar _tmp31 = std::asinh(_tmp30);
  const Scalar _tmp32 = Scalar(9.6622558468725703) * _tmp31;
  const Scalar _tmp33 =
      -Scalar(0.1034955) * _tmp32 * fh1 -
      Scalar(0.50022801989500498) *
          std::sqrt(
              Scalar(std::pow(Scalar(1 - Scalar(0.20689664689659551) * _tmp29), Scalar(2)) +
                     Scalar(0.13817235445745474) *
                         std::pow(Scalar(-Scalar(0.55659957866191134) * _tmp19 - 1), Scalar(2))));
  const Scalar _tmp34 = _tmp0 * _tmp33;
  const Scalar _tmp35 = std::pow(fh1, Scalar(-2));
  const Scalar _tmp36 =
      std::pow(Scalar(_tmp35 * std::pow(fv1, Scalar(2)) + 1), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp37 = Scalar(1.0) * _tmp31;
  const Scalar _tmp38 = _tmp29 + Scalar(-4.8333311099999996);
  const Scalar _tmp39 = _tmp19 + Scalar(1.79662371);
  const Scalar _tmp40 = std::pow(Scalar(std::pow(_tmp38, Scalar(2)) + std::pow(_tmp39, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp41 = _tmp38 * _tmp40;
  const Scalar _tmp42 = -Scalar(0.20999999999999999) * _tmp10 + Scalar(0.20999999999999999) * _tmp8;
  const Scalar _tmp43 = -_tmp42;
  const Scalar _tmp44 = -Scalar(0.010999999999999999) * _tmp2 -
                        Scalar(0.010999999999999999) * _tmp20 + Scalar(-0.010999999999999999);
  const Scalar _tmp45 = Scalar(0.20999999999999999) * _tmp23 + Scalar(0.20999999999999999) * _tmp24;
  const Scalar _tmp46 = _tmp44 + _tmp45;
  const Scalar _tmp47 = _tmp43 + _tmp46;
  const Scalar _tmp48 = _tmp47 * fh1;
  const Scalar _tmp49 = -Scalar(3.29616) * _tmp25 - _tmp28 * fv1 - _tmp41 * _tmp48;
  const Scalar _tmp50 = _tmp12 - _tmp16;
  const Scalar _tmp51 = _tmp5 + _tmp50;
  const Scalar _tmp52 = _tmp50 + _tmp6;
  const Scalar _tmp53 = Scalar(1.0) * _tmp52;
  const Scalar _tmp54 = -_tmp21;
  const Scalar _tmp55 = _tmp27 + _tmp54;
  const Scalar _tmp56 = Scalar(1.0) * _tmp55;
  const Scalar _tmp57 = -_tmp56;
  const Scalar _tmp58 = _tmp22 + _tmp26;
  const Scalar _tmp59 = _tmp54 + _tmp58;
  const Scalar _tmp60 = _tmp57 + _tmp59;
  const Scalar _tmp61 = _tmp21 + _tmp58;
  const Scalar _tmp62 = Scalar(1.0) / (_tmp57 + _tmp61);
  const Scalar _tmp63 = _tmp17 + _tmp5;
  const Scalar _tmp64 = _tmp62 * (_tmp53 - _tmp63);
  const Scalar _tmp65 = _tmp60 * _tmp64;
  const Scalar _tmp66 = Scalar(1.0) / (-_tmp51 + _tmp53 - _tmp65);
  const Scalar _tmp67 = Scalar(1.0) * _tmp66;
  const Scalar _tmp68 = _tmp64 * _tmp67;
  const Scalar _tmp69 = _tmp62 * (_tmp65 * _tmp67 + Scalar(1.0));
  const Scalar _tmp70 = Scalar(1.0) * _tmp68 - Scalar(1.0) * _tmp69;
  const Scalar _tmp71 = _tmp53 + _tmp56 * _tmp64;
  const Scalar _tmp72 = _tmp66 * _tmp71;
  const Scalar _tmp73 = _tmp44 - _tmp45;
  const Scalar _tmp74 = _tmp42 + _tmp73;
  const Scalar _tmp75 = _tmp51 + p_init0;
  const Scalar _tmp76 = _tmp75 + Scalar(-2.5202214700000001);
  const Scalar _tmp77 = _tmp59 + p_init1;
  const Scalar _tmp78 = _tmp77 + Scalar(8.3888750099999996);
  const Scalar _tmp79 = std::pow(Scalar(std::pow(_tmp76, Scalar(2)) + std::pow(_tmp78, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp80 = _tmp78 * _tmp79;
  const Scalar _tmp81 = _tmp43 + _tmp73;
  const Scalar _tmp82 = _tmp55 + p_init1;
  const Scalar _tmp83 = _tmp82 + Scalar(8.3196563700000006);
  const Scalar _tmp84 = _tmp52 + p_init0;
  const Scalar _tmp85 = _tmp84 + Scalar(1.9874742000000001);
  const Scalar _tmp86 = Scalar(1.0) / (_tmp85);
  const Scalar _tmp87 = _tmp83 * _tmp86;
  const Scalar _tmp88 = _tmp81 * _tmp87;
  const Scalar _tmp89 = _tmp63 + p_init0;
  const Scalar _tmp90 = _tmp89 + Scalar(-2.71799795);
  const Scalar _tmp91 = _tmp61 + p_init1;
  const Scalar _tmp92 = _tmp91 + Scalar(-4.7752063900000001);
  const Scalar _tmp93 = std::pow(Scalar(std::pow(_tmp90, Scalar(2)) + std::pow(_tmp92, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp94 = _tmp90 * _tmp93;
  const Scalar _tmp95 = _tmp42 + _tmp46;
  const Scalar _tmp96 = _tmp92 * _tmp93;
  const Scalar _tmp97 = -_tmp88 * _tmp94 + _tmp95 * _tmp96;
  const Scalar _tmp98 = _tmp76 * _tmp79;
  const Scalar _tmp99 = -_tmp80 + _tmp87 * _tmp98;
  const Scalar _tmp100 = Scalar(1.0) / (_tmp87 * _tmp94 - _tmp96);
  const Scalar _tmp101 = _tmp100 * _tmp99;
  const Scalar _tmp102 = _tmp81 * _tmp98;
  const Scalar _tmp103 = -_tmp101 * _tmp97 - _tmp102 * _tmp87 + _tmp74 * _tmp80;
  const Scalar _tmp104 = _tmp81 * _tmp94 - _tmp94 * _tmp95;
  const Scalar _tmp105 = -_tmp101 * _tmp104 + _tmp102 - _tmp103 * _tmp64 - _tmp74 * _tmp98;
  const Scalar _tmp106 = Scalar(1.0) / (_tmp105);
  const Scalar _tmp107 = 0;
  const Scalar _tmp108 = -_tmp103 * _tmp107 + _tmp57 - _tmp60 * _tmp72;
  const Scalar _tmp109 = Scalar(1.0) * _tmp62;
  const Scalar _tmp110 = Scalar(43.164000000000001) - fv1;
  const Scalar _tmp111 = Scalar(1.0) * _tmp100;
  const Scalar _tmp112 = _tmp111 * _tmp97;
  const Scalar _tmp113 = -_tmp104 * _tmp111 + _tmp112 * _tmp64;
  const Scalar _tmp114 =
      std::sqrt(Scalar(std::pow(_tmp83, Scalar(2)) + std::pow(_tmp85, Scalar(2))));
  const Scalar _tmp115 = Scalar(1.0) / (_tmp114);
  const Scalar _tmp116 = _tmp114 * _tmp86;
  const Scalar _tmp117 = _tmp116 * (_tmp115 * _tmp52 * _tmp83 - _tmp115 * _tmp55 * _tmp85);
  const Scalar _tmp118 = _tmp100 * (_tmp117 * _tmp94 + _tmp61 * _tmp94 - _tmp63 * _tmp96);
  const Scalar _tmp119 = _tmp117 * _tmp98 - _tmp118 * _tmp99 - _tmp51 * _tmp80 + _tmp59 * _tmp98;
  const Scalar _tmp120 = _tmp106 * _tmp119;
  const Scalar _tmp121 = Scalar(1.0) / (_tmp119);
  const Scalar _tmp122 = _tmp105 * _tmp121;
  const Scalar _tmp123 = _tmp122 * (-_tmp113 * _tmp120 - Scalar(1.0) * _tmp118);
  const Scalar _tmp124 = _tmp60 * _tmp66;
  const Scalar _tmp125 = _tmp106 * (_tmp113 + _tmp123);
  const Scalar _tmp126 = -_tmp103 * _tmp125 - _tmp112 + _tmp123 * _tmp124;
  const Scalar _tmp127 = Scalar(1.0) * _tmp41 * (-_tmp109 * _tmp126 + _tmp123 * _tmp67);
  const Scalar _tmp128 = _tmp100 * _tmp87;
  const Scalar _tmp129 = _tmp128 * _tmp97 + _tmp88;
  const Scalar _tmp130 = _tmp104 * _tmp128 - _tmp129 * _tmp64 - _tmp81;
  const Scalar _tmp131 = _tmp122 * (-_tmp117 + _tmp118 * _tmp87 - _tmp120 * _tmp130);
  const Scalar _tmp132 = _tmp106 * (_tmp130 + _tmp131);
  const Scalar _tmp133 = -_tmp103 * _tmp132 + _tmp124 * _tmp131 + _tmp129;
  const Scalar _tmp134 = _tmp39 * _tmp40;
  const Scalar _tmp135 = Scalar(1.0) * _tmp134 * (-_tmp109 * _tmp133 + _tmp131 * _tmp67);
  const Scalar _tmp136 = _tmp122 * _tmp67;
  const Scalar _tmp137 = Scalar(1.0) * _tmp121;
  const Scalar _tmp138 = -_tmp103 * _tmp137 + _tmp136 * _tmp60;
  const Scalar _tmp139 = -_tmp134 * _tmp28 + _tmp18 * _tmp41;
  const Scalar _tmp140 = Scalar(1.0) * _tmp139 * (-_tmp109 * _tmp138 + _tmp136);
  const Scalar _tmp141 = Scalar(3.29616) * _tmp11 + _tmp134 * _tmp48 + _tmp18 * fv1;
  const Scalar _tmp142 = _tmp60 * _tmp62;
  const Scalar _tmp143 = _tmp142 * _tmp67;
  const Scalar _tmp144 = Scalar(1.0) * _tmp143 - Scalar(1.0) * _tmp67;
  const Scalar _tmp145 =
      Scalar(1.0) * _tmp110 * (-_tmp108 * _tmp109 - _tmp67 * _tmp71 + Scalar(1.0)) + _tmp127 * fh1 +
      _tmp135 * fh1 + _tmp140 * fh1 + _tmp141 * _tmp144 + _tmp49 * _tmp70;
  const Scalar _tmp146 = -_tmp125 * _tmp99 + Scalar(1.0);
  const Scalar _tmp147 = _tmp100 * _tmp94;
  const Scalar _tmp148 = _tmp116 * _tmp41 * (_tmp125 * _tmp98 + _tmp146 * _tmp147);
  const Scalar _tmp149 = -_tmp132 * _tmp99 - _tmp87;
  const Scalar _tmp150 = _tmp116 * _tmp134 * (_tmp132 * _tmp98 + _tmp147 * _tmp149 + Scalar(1.0));
  const Scalar _tmp151 = _tmp111 * _tmp121 * _tmp99;
  const Scalar _tmp152 = _tmp116 * _tmp139 * (_tmp137 * _tmp98 - _tmp151 * _tmp94);
  const Scalar _tmp153 = -_tmp110 * _tmp116 * (-_tmp101 * _tmp107 * _tmp94 + _tmp107 * _tmp98) -
                         _tmp148 * fh1 - _tmp150 * fh1 - _tmp152 * fh1;
  const Scalar _tmp154 = Scalar(1.0) / (_tmp153);
  const Scalar _tmp155 = std::asinh(_tmp145 * _tmp154);
  const Scalar _tmp156 = Scalar(9.6622558468725703) * _tmp153;
  const Scalar _tmp157 =
      -_tmp155 * _tmp156 -
      Scalar(8.3196563700000006) *
          std::sqrt(
              Scalar(std::pow(Scalar(-Scalar(0.12019727204189803) * _tmp82 - 1), Scalar(2)) +
                     Scalar(0.057067943376852184) *
                         std::pow(Scalar(-Scalar(0.50315118556004401) * _tmp84 - 1), Scalar(2))));
  const Scalar _tmp158 = Scalar(0.1034955) * _tmp154;
  const Scalar _tmp159 = _tmp157 * _tmp158;
  const Scalar _tmp160 = Scalar(1.0) * _tmp155;
  const Scalar _tmp161 = -_tmp148 - _tmp150 - _tmp152;
  const Scalar _tmp162 = Scalar(9.6622558468725703) * _tmp161;
  const Scalar _tmp163 = std::pow(_tmp153, Scalar(-2));
  const Scalar _tmp164 = _tmp41 * _tmp47;
  const Scalar _tmp165 = _tmp134 * _tmp47;
  const Scalar _tmp166 = _tmp161 * _tmp163;
  const Scalar _tmp167 = (-_tmp145 * _tmp166 + _tmp154 * (_tmp127 + _tmp135 + _tmp140 +
                                                          _tmp144 * _tmp165 - _tmp164 * _tmp70)) /
                         std::sqrt(Scalar(std::pow(_tmp145, Scalar(2)) * _tmp163 + 1));
  const Scalar _tmp168 = _tmp133 * _tmp134 * _tmp62;
  const Scalar _tmp169 = _tmp126 * _tmp41 * _tmp62;
  const Scalar _tmp170 = _tmp141 * _tmp67;
  const Scalar _tmp171 = _tmp138 * _tmp139 * _tmp62;
  const Scalar _tmp172 = _tmp108 * _tmp110 * _tmp62 - _tmp142 * _tmp170 + _tmp168 * fh1 +
                         _tmp169 * fh1 + _tmp171 * fh1 + _tmp49 * _tmp69;
  const Scalar _tmp173 = _tmp107 * _tmp110;
  const Scalar _tmp174 = _tmp100 * _tmp134 * _tmp149;
  const Scalar _tmp175 = _tmp100 * _tmp146 * _tmp41;
  const Scalar _tmp176 = _tmp139 * _tmp151;
  const Scalar _tmp177 = -_tmp101 * _tmp173 + _tmp174 * fh1 + _tmp175 * fh1 - _tmp176 * fh1;
  const Scalar _tmp178 = Scalar(1.0) / (_tmp177);
  const Scalar _tmp179 = std::asinh(_tmp172 * _tmp178);
  const Scalar _tmp180 = Scalar(1.0) * _tmp179;
  const Scalar _tmp181 = Scalar(9.6622558468725703) * _tmp177;
  const Scalar _tmp182 =
      -_tmp179 * _tmp181 -
      Scalar(4.7752063900000001) *
          std::sqrt(
              Scalar(Scalar(0.32397683292140877) *
                         std::pow(Scalar(1 - Scalar(0.36791786395571047) * _tmp89), Scalar(2)) +
                     std::pow(Scalar(1 - Scalar(0.20941503221602112) * _tmp91), Scalar(2))));
  const Scalar _tmp183 = Scalar(0.1034955) * _tmp178;
  const Scalar _tmp184 = _tmp182 * _tmp183;
  const Scalar _tmp185 = _tmp174 + _tmp175 - _tmp176;
  const Scalar _tmp186 = Scalar(9.6622558468725703) * _tmp185;
  const Scalar _tmp187 = std::pow(_tmp177, Scalar(-2));
  const Scalar _tmp188 = _tmp185 * _tmp187;
  const Scalar _tmp189 = (-_tmp172 * _tmp188 + _tmp178 * (-_tmp143 * _tmp165 - _tmp164 * _tmp69 +
                                                          _tmp168 + _tmp169 + _tmp171)) /
                         std::sqrt(Scalar(std::pow(_tmp172, Scalar(2)) * _tmp187 + 1));
  const Scalar _tmp190 = _tmp136 * _tmp139;
  const Scalar _tmp191 = _tmp123 * _tmp41 * _tmp66;
  const Scalar _tmp192 = _tmp131 * _tmp134 * _tmp66;
  const Scalar _tmp193 =
      _tmp110 * _tmp72 + _tmp170 - _tmp190 * fh1 - _tmp191 * fh1 - _tmp192 * fh1 - _tmp49 * _tmp68;
  const Scalar _tmp194 = _tmp125 * _tmp41;
  const Scalar _tmp195 = _tmp137 * _tmp139;
  const Scalar _tmp196 = _tmp132 * _tmp134;
  const Scalar _tmp197 = _tmp173 + _tmp194 * fh1 + _tmp195 * fh1 + _tmp196 * fh1;
  const Scalar _tmp198 = Scalar(1.0) / (_tmp197);
  const Scalar _tmp199 = std::asinh(_tmp193 * _tmp198);
  const Scalar _tmp200 = Scalar(1.0) * _tmp199;
  const Scalar _tmp201 = std::pow(_tmp197, Scalar(-2));
  const Scalar _tmp202 = _tmp194 + _tmp195 + _tmp196;
  const Scalar _tmp203 = _tmp201 * _tmp202;
  const Scalar _tmp204 = (-_tmp193 * _tmp203 + _tmp198 * (_tmp164 * _tmp68 + _tmp165 * _tmp67 -
                                                          _tmp190 - _tmp191 - _tmp192)) /
                         std::sqrt(Scalar(std::pow(_tmp193, Scalar(2)) * _tmp201 + 1));
  const Scalar _tmp205 = Scalar(9.6622558468725703) * _tmp197;
  const Scalar _tmp206 =
      -_tmp199 * _tmp205 -
      Scalar(8.3888750099999996) *
          std::sqrt(
              Scalar(Scalar(0.090254729040973036) *
                         std::pow(Scalar(1 - Scalar(0.39679052492160538) * _tmp75), Scalar(2)) +
                     std::pow(Scalar(-Scalar(0.11920549523123722) * _tmp77 - 1), Scalar(2))));
  const Scalar _tmp207 = Scalar(0.1034955) * _tmp198;
  const Scalar _tmp208 = _tmp206 * _tmp207;
  const Scalar _tmp209 = Scalar(9.6622558468725703) * _tmp202;

  // Output terms (1)
  Eigen::Matrix<Scalar, 4, 1> _res;

  _res(0, 0) =
      Scalar(9.6622558468725703) * fh1 *
          (Scalar(1.0) * _tmp35 * _tmp36 * fv1 * std::cosh(_tmp37) -
           (Scalar(0.1034955) * _tmp0 * (Scalar(9.6622558468725703) * _tmp30 * _tmp36 - _tmp32) -
            _tmp33 * _tmp35) *
               std::cosh(_tmp34)) -
      Scalar(9.6622558468725703) * std::sinh(_tmp34) -
      Scalar(9.6622558468725703) * std::sinh(_tmp37);
  _res(1, 0) = _tmp156 * (-Scalar(1.0) * _tmp167 * std::cosh(_tmp160) -
                          (-Scalar(0.1034955) * _tmp157 * _tmp166 +
                           _tmp158 * (-_tmp155 * _tmp162 - _tmp156 * _tmp167)) *
                              std::cosh(_tmp159)) +
               _tmp162 * (-std::sinh(_tmp159) - std::sinh(_tmp160));
  _res(2, 0) = _tmp181 * (-Scalar(1.0) * _tmp189 * std::cosh(_tmp180) -
                          (-Scalar(0.1034955) * _tmp182 * _tmp188 +
                           _tmp183 * (-_tmp179 * _tmp186 - _tmp181 * _tmp189)) *
                              std::cosh(_tmp184)) +
               _tmp186 * (-std::sinh(_tmp180) - std::sinh(_tmp184));
  _res(3, 0) = _tmp205 * (-Scalar(1.0) * _tmp204 * std::cosh(_tmp200) -
                          (-Scalar(0.1034955) * _tmp203 * _tmp206 +
                           _tmp207 * (-_tmp199 * _tmp209 - _tmp204 * _tmp205)) *
                              std::cosh(_tmp208)) +
               _tmp209 * (-std::sinh(_tmp200) - std::sinh(_tmp208));

  return _res;
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym

# Paper supplement

Every PDF here carries a caption that explains what the figure or table shows and how to read it.

## Figures (`figures/`)

Trail figures on the subround diagrams. Notation as in Section 3.6 of the paper: red italic bits are always active, green bits vary over the kept trails, `[ ]` means no mask.

- `5round-f1-bit-propagation.pdf`: Trail of f_1 in the 5-round attack: subrounds 19 and 20. Notation as in Section 3.6 of the paper.
- `5.5round-f1-bit-propagation.pdf`: Trail of f_1 in the 5.5-round attack: subrounds 19 to 22. Notation as in Section 3.6 of the paper.
- `5.5round-f1-final-key-addition.pdf`: Trail of f_1 in the 5.5-round attack: final key addition bar z_i = v_i^{[22]} boxplus k_i. Notation as in Section 3.6 of the paper.
- `5.5round-f2-bit-propagation.pdf`: Trail of f_2 in the 5.5-round attack: subround 22. Notation as in Section 3.6 of the paper.
- `5.5round-f2-final-key-addition.pdf`: Trail of f_2 in the 5.5-round attack: final key addition bar z_i = v_i^{[22]} boxplus k_i. Notation as in Section 3.6 of the paper.
- `5round-mask-extension-subrounds-17-20.pdf`: Extension of the output mask through subrounds 17 to 20 (5-round attack). Blue: a bit of the final parity. Black: an intermediate bit that a later subround uses. A masked word that the next subrounds do not use keeps its value, so it is labeled with the index it has when it is next used, or 20 if it is not used again. Notation as in Section 3.6 of the paper.
- `5.5round-mask-extension-subrounds-17-22.pdf`: Extension of the output mask through subrounds 17 to 22 (5.5-round attack). Top row: subrounds 17 to 19. Bottom row: subrounds 20 to 22. Blue: a bit of the final parity. Black: an intermediate bit that a later subround uses. A masked word that the next subrounds do not use keeps its value, so it is labeled with the index it has when it is next used, or 22 if it is not used again. Notation as in Section 3.6 of the paper.
- `6round-mask-extension-subrounds-17-24.pdf`: Extension of the output mask through subrounds 17 to 24 (6-round attack). The rows show subrounds 17 to 19, 20 to 22, and 23 and 24. Blue: a bit of the final parity. Black: an intermediate bit that a later subround uses. A masked word that the next subrounds do not use keeps its value, so it is labeled with the index it has when it is next used, or 24 if it is not used again. Notation as in Section 3.6 of the paper.
- `6.5round-key-bit-classification.pdf`: Key bits of the 6.5-round attack, each word shown as k_i[31-0] from left to right. Gray: the 113 CPNBs. Red: the restricted bits k^res. Blue: the syncopated segments. White: the other significant key bits. The removed starred condition (k_9[7] and its segment) is still shown.
- `6.75round-key-bit-classification.pdf`: Key bits of the 6.75-round attack, each word shown as k_i[31-0] from left to right. Gray: the 80 CPNBs. Red: the restricted bits k^res. Blue: the syncopated segments. White: the other significant key bits. The two removed starred conditions (k_2[14] and k_9[8], with their segments) are still shown.

## Tables (`tables/`)

`tables.pdf` holds six tables:

- the per-component PNB sets P_i of each single-bit component of the output mask, with their bit indices, for the 5-, 5.5- and 6-round attacks (Appendix B of the paper);
- the full trail enumerations of f_1, f_2 and f_3 in the 5.5-round attack.

Key bits are numbered as in Equation 1 of the paper.

import { micromark } from "micromark";
import {
  gfmAutolinkLiteral,
  gfmAutolinkLiteralHtml,
} from "micromark-extension-gfm-autolink-literal";

export const parse = (content = "") =>
  micromark(content, {
    extensions: [gfmAutolinkLiteral()],
    htmlExtensions: [gfmAutolinkLiteralHtml()],
  })
    .replaceAll("<p>", "")
    .replaceAll("</p>", "");

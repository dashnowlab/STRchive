/** @type {import("prettier").Config} */
export default {
  plugins: [
    "prettier-plugin-astro",
    "@ianvs/prettier-plugin-sort-imports",
    "prettier-plugin-css-order",
  ],
  importOrder: [
    "^react",
    "^[a-zA-Z]",
    "^@[a-zA-Z]",
    "^@/",
    "^/",
    "^./",
    "^../",
  ],
  importOrderParserPlugins: ["jsx", "importAssertions"],
  cssDeclarationSorterOrder: "smacss",
  overrides: [
    {
      files: "*.astro",
      options: {
        parser: "astro",
      },
    },
  ],
};

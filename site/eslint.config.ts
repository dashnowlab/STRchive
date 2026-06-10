import js from "@eslint/js";
import astro from "eslint-plugin-astro";
import tailwind from "eslint-plugin-better-tailwindcss";
import * as mdx from "eslint-plugin-mdx";
import prettier from "eslint-plugin-prettier/recommended";
import reactHooks from "eslint-plugin-react-hooks";
import { defineConfig, globalIgnores } from "eslint/config";
import globals from "globals";
import tslint from "typescript-eslint";

export default defineConfig([
  globalIgnores(["dist", "public"]),

  // https://github.com/mdx-js/eslint-mdx/issues/92
  {
    name: "TypeScript",
    extends: tslint.configs.recommended,
    ignores: ["**/*.mdx"],
    rules: {
      "@typescript-eslint/no-unused-vars": ["warn", { caughtErrors: "none" }],
      "@typescript-eslint/consistent-type-definitions": ["error", "type"],
      "@typescript-eslint/consistent-type-imports": "error",
    },
    languageOptions: {
      parserOptions: {
        project: true,
        tsconfigRootDir: import.meta.dirname,
      },
    },
  },

  {
    name: "Astro",
    extends: astro.configs.recommended,
  },

  {
    name: "JavaScript",
    files: ["**/*.{js,jsx,cjs,mjs}"],
    ...js.configs.recommended,
    rules: {
      "prefer-const": ["error", { destructuring: "all" }],
    },
  },

  {
    name: "React Hooks",
    extends: [reactHooks.configs.flat.recommended],
  },

  {
    name: "Prettier",
    extends: [prettier],
    ignores: ["**/*.mdx"],
    rules: {
      "prettier/prettier": "warn",
    },
  },

  {
    name: "Tailwind",
    files: ["**/*.{js,jsx,ts,tsx}"],
    extends: [tailwind.configs.recommended],
    rules: {
      "better-tailwindcss/enforce-consistent-line-wrapping": [
        "warn",
        {
          preferSingleLine: true,
          group: "never",
          printWidth: 0,
        },
      ],
      "better-tailwindcss/no-unknown-classes": ["warn"],
    },
    settings: {
      "better-tailwindcss": { entryPoint: "./src/styles.css" },
    },
  },

  {
    name: "MDX",
    ...mdx.flat,
    rules: {
      "@typescript-eslint/consistent-type-imports": "off",
    },
  },

  {
    languageOptions: {
      globals: globals.browser,
      ecmaVersion: 2020,
    },
  },
]);

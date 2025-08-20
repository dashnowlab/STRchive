import Ajv from "ajv/dist/2020";
import schema from "~/STRchive-loci.schema.json";

for (const [key, { type }] of Object.entries(schema.properties))
  if ([type].flat().includes("null")) {
    const index = schema.required.indexOf(key);
    if (index !== -1) schema.required.splice(index, 1);
  }

export const ajv = new Ajv({ strict: false, allErrors: true });
export const validator = ajv.compile(schema);
export const validate = (data) => {
  validator(data);
  return validator.errors;
};

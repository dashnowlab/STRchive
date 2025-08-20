import Ajv from "ajv/dist/2020";
import schema from "~/STRchive-loci.schema.json";

export const ajv = new Ajv({ strict: false, allErrors: true });
export const validator = ajv.compile(schema);
export const validate = (data) => {
  validator(data);
  return validator.errors;
};
